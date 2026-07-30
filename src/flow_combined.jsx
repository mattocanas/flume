import { useState, useRef, useCallback, useMemo, useEffect } from "react";
import { parseFCS, detectScatter } from "./fcs.js";

// ─── Safe min/max ────────────────────────────────────
function arrMin(a){let m=Infinity;for(let i=0;i<a.length;i++)if(a[i]<m)m=a[i];return m;}
function arrMax(a){let m=-Infinity;for(let i=0;i<a.length;i++)if(a[i]>m)m=a[i];return m;}

// ─── Transform (asinh ≈ biexponential) ───────────────
const COF=150;
const T=x=>Math.asinh(x/COF);
const invT=y=>Math.sinh(y)*COF;
const N_BINS=220;
const SUP={0:"⁰",1:"¹",2:"²",3:"³",4:"⁴",5:"⁵",6:"⁶",7:"⁷",8:"⁸",9:"⁹"};

// Five tasteful, mutually distinct qualitative palettes. Each is ordered so adjacent samples
// (which take adjacent colors) stay far apart in color — validated with the dataviz palette
// checker (normal-vision ΔE floor). "Colorblind" additionally clears the CVD separation gate.
const PALETTES={
  Classic:["#4E79A7","#F28E2B","#E15759","#76B7B2","#EDC948","#59A14F","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"],
  Bright:["#E58606","#5D69B1","#99C945","#CC61B0","#52BCA3","#DAA51B","#2F8AC4","#764E9F","#ED645A","#CC3A8E"],
  Bold:["#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72"],
  Spectrum:["#2E6FB0","#E4572E","#17A398","#F2B705","#8E44AD","#3FA34D","#E5679E","#8C6D31","#2CBDCE","#D1495B","#6A8D2F","#5C5CD6"],
  Colorblind:["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#111111"],
};
const PALETTE=PALETTES.Classic;
const COLOR_OPTIONS=[...new Set(Object.values(PALETTES).flat())];

const PL={w:420,h:290,ml:58,mr:22,mt:34,mb:50};
const PW=PL.w-PL.ml-PL.mr;
const PH=PL.h-PL.mt-PL.mb;

// ─── Utilities ───────────────────────────────────────
function parseCSV(text){
  const lines=text.trim().split("\n").filter(l=>l.trim());
  const sep=lines[0].includes("\t")?"\t":",";
  const headers=lines[0].split(sep).map(h=>h.trim().replace(/^"|"$/g,""));
  const columns={};headers.forEach(h=>(columns[h]=[]));
  for(let i=1;i<lines.length;i++){const vals=lines[i].split(sep);headers.forEach((h,j)=>{const v=parseFloat(vals[j]);if(!isNaN(v))columns[h].push(v);});}
  return{headers,columns};
}
function detectPE(headers){
  const pats=["PE-A","YL1-A","R-PE-A","PE-H","YL1-H","BL2-A","PE-Cy5-A","PE-Cy7-A"];
  for(const p of pats)if(headers.includes(p))return p;
  const f=headers.find(h=>/\bPE\b/i.test(h));return f||headers.find(h=>!/FSC|SSC|Time/i.test(h))||headers[0];
}
function detectChannels(headers){
  const fluoro=headers.filter(h=>!/FSC|SSC|Time/i.test(h));
  return fluoro.length>=2?[fluoro[0],fluoro[1]]:[headers[0],headers[1]||headers[0]];
}
function computeHist(values,nBins,lo,hi){
  const w=(hi-lo)/nBins;const counts=new Array(nBins).fill(0);
  const centers=Array.from({length:nBins},(_,i)=>lo+(i+0.5)*w);
  for(let k=0;k<values.length;k++){const tv=T(values[k]);const idx=Math.floor((tv-lo)/w);if(idx>=0&&idx<nBins)counts[idx]++;}
  return{counts,centers};
}
function fmtTick(v){
  if(v===0)return "0";
  const s=v<0?"−":"";const a=Math.abs(v);
  if(a<1)return "0";
  if(a<100)return s+Math.round(a);
  const e=Math.floor(Math.log10(a));
  const m=a/Math.pow(10,e);
  const sup=SUP[e]||"⁸⁻"+SUP[e%10];
  const mRound=Math.round(m*10)/10;
  if(Math.abs(mRound-1)<0.05)return s+"10"+sup;
  return s+mRound+"×10"+sup;
}
// Pixel-aware biexponential ticks. Returns [{v,label}]: marks thinned so they never
// bunch up near 0; only "major" values (powers of 10 and 5×10ⁿ) get labels, as space allows.
function axisTicks(dMin,dMax,pxW){
  const tL=T(dMin),tH=T(dMax);const span=(tH-tL)||1;const p=span*0.03;
  const cand=new Set([0]);
  for(let e=0;e<=6;e++){const pow=10**e;for(const m of[1,2,5])cand.add(m*pow);} // positive ticks only
  const inRange=[...cand].sort((a,b)=>a-b).filter(x=>T(x)>=tL-p&&T(x)<=tH+p);
  const px=v=>(T(v)-tL)/span*pxW;
  const isMajor=v=>{if(v===0)return true;const a=Math.abs(v);const l=Math.log10(a);const fl=Math.floor(l+1e-9);return Math.abs(l-Math.round(l))<1e-6||Math.abs(a/(10**fl)-5)<1e-6;};
  const ticks=[];let lastX=-1e9;
  for(const v of inRange){const x=px(v);if(v===0||x-lastX>=15){ticks.push(v);lastX=x;}}
  const labels=new Set();let lastL=-1e9;
  for(const v of ticks){if(!isMajor(v))continue;const x=px(v);if(v===0||x-lastL>=30){labels.add(v);lastL=x;}}
  return ticks.map(v=>({v,label:labels.has(v)}));
}
// Classic log axis: decade labels (10ⁿ) + log-spaced minor ticks (2–9 per decade). Positions in log10.
function logTicks(dMin,dMax,pxW){
  const lL=Math.log10(dMin),lH=Math.log10(dMax);const span=(lH-lL)||1;const p=span*0.02;
  const px=v=>(Math.log10(v)-lL)/span*pxW;
  const cand=[];
  for(let e=Math.floor(lL)-1;e<=Math.ceil(lH)+1;e++)for(let m=1;m<=9;m++){const v=m*10**e;const lv=Math.log10(v);if(lv>=lL-p&&lv<=lH+p)cand.push(v);}
  cand.sort((a,b)=>a-b);
  const isPow=v=>{const e=Math.log10(v);return Math.abs(e-Math.round(e))<1e-9;};
  const ticks=[];let lastX=-1e9;
  for(const v of cand){const x=px(v);if(isPow(v)||x-lastX>=4.5){ticks.push(v);lastX=x;}}
  // When the axis spans only 1–2 decades, decade powers alone give too few labels, so also
  // label the 2× and 5× sub-decade ticks (e.g. 2×10⁵, 5×10⁵) as on classic FACS axes.
  const nPow=ticks.filter(isPow).length;
  const subOK=nPow<=3;
  const isLabelCand=v=>{if(isPow(v))return true;if(!subOK)return false;const e=Math.floor(Math.log10(v)+1e-9);const m=Math.round(v/10**e);return m===2||m===5;};
  const labels=new Set();let lastL=-1e9;
  for(const v of ticks){if(!isLabelCand(v))continue;const x=px(v);if(x-lastL>=18){labels.add(v);lastL=x;}}
  return ticks.map(v=>({v,label:labels.has(v)}));
}
function fmtLog(v){const e=Math.round(Math.log10(v));return "10"+(SUP[e]||("^"+e));}
// Log tick label that also renders sub-decade values as m×10ⁿ (for 1–2 decade axes)
function fmtLogTick(v){const e=Math.floor(Math.log10(v)+1e-9);const m=Math.round(v/10**e);if(m>=10)return "10"+(SUP[e+1]||("^"+(e+1)));if(m<=1)return "10"+(SUP[e]||("^"+e));return m+"×10"+(SUP[e]||("^"+e));}
// Separable box blur of a density grid — turns blocky per-bin counts into a smooth KDE-like field
// so cluster colors transition as clean gradients instead of speckle.
function boxBlur(src,w,h,r,passes){
  let a=Float32Array.from(src);
  for(let p=0;p<passes;p++){
    const b=new Float32Array(a.length);
    for(let y=0;y<h;y++)for(let x=0;x<w;x++){let s=0,c=0;for(let dx=-r;dx<=r;dx++){const xx=x+dx;if(xx<0||xx>=w)continue;s+=a[y*w+xx];c++;}b[y*w+x]=s/c;}
    const a2=new Float32Array(a.length);
    for(let y=0;y<h;y++)for(let x=0;x<w;x++){let s=0,c=0;for(let dy=-r;dy<=r;dy++){const yy=y+dy;if(yy<0||yy>=h)continue;s+=b[yy*w+x];c++;}a2[y*w+x]=s/c;}
    a=a2;
  }
  return a;
}
// Paint density-colored dots into ctx over the plot rect (x0,y0,PW,PH). Smooths the density field,
// then draws events sorted by density so the dense (warm) points sit on TOP of the sparse blue halo.
function paintDensity(ctx,xVals,yVals,xS,yS,x0,y0,PW,PH,dotSize){
  const n=Math.min(xVals.length,yVals.length);
  const BS=2,gc=Math.ceil(PW/BS),gr=Math.ceil(PH/BS);
  const grid=new Float32Array(gc*gr),bx=new Int32Array(n);
  for(let i=0;i<n;i++){const gx=(xS(xVals[i])-x0)/BS|0,gy=(yS(yVals[i])-y0)/BS|0;if(gx<0||gy<0||gx>=gc||gy>=gr){bx[i]=-1;}else{const b=gy*gc+gx;bx[i]=b;grid[b]++;}}
  const sm=boxBlur(grid,gc,gr,2,2);
  const nz=[];for(let i=0;i<sm.length;i++)if(sm[i]>0)nz.push(sm[i]);nz.sort((a,b)=>a-b);
  const refD=nz.length?nz[Math.min(nz.length-1,Math.floor(nz.length*0.995))]:1,lr=Math.log(1+refD)||1;
  const step=n>60000?Math.ceil(n/60000):1;
  const idx=[];for(let i=0;i<n;i+=step)if(bx[i]>=0)idx.push(i);
  idx.sort((a,b)=>sm[bx[a]]-sm[bx[b]]); // draw sparse first, dense last (painter's algorithm)
  const ds=Math.max(1,Math.round(dotSize)),off=(ds-1)>>1;
  for(let k=0;k<idx.length;k++){const i=idx[k];const t=Math.min(1,Math.log(1+sm[bx[i]])/lr);ctx.fillStyle=DENSITY[Math.min(63,Math.floor(Math.pow(t,0.9)*63))];ctx.fillRect((xS(xVals[i])|0)-off,(yS(yVals[i])|0)-off,ds,ds);}
}
// Auto log range snapped to decade bounds
function logRange(values){
  const pos=values.filter(v=>v>0).sort((a,b)=>a-b);
  if(pos.length<3)return{lo:0,hi:5,dMin:1,dMax:100000};
  const lo=pos[Math.floor(pos.length*0.003)];
  const hi=pos[Math.min(pos.length-1,Math.floor(pos.length*0.9997))];
  const eLo=Math.max(0,Math.floor(Math.log10(lo)));
  const eHi=Math.max(eLo+1,Math.ceil(Math.log10(hi)));
  const dMin=10**eLo,dMax=10**eHi;
  return{lo:Math.log10(dMin),hi:Math.log10(dMax),dMin,dMax};
}
function manualLog(r,auto){
  if(!auto)return auto;const mn=Number(r.min),mx=Number(r.max);
  if(!r.min.trim()||!r.max.trim()||!isFinite(mn)||!isFinite(mx)||mx<=mn||mn<=0)return auto;
  return{lo:Math.log10(mn),hi:Math.log10(mx),dMin:mn,dMax:mx};
}
function getTicks(dMin,dMax){
  const tL=T(dMin),tH=T(dMax);const span=tH-tL;const p=span*0.03;
  const candidates=new Set([0]);
  for(let e=0;e<=6;e++){
    const pow=10**e;
    for(const m of [1,2,5]){
      candidates.add(m*pow);
      candidates.add(-m*pow);
    }
  }
  let inRange=[...candidates].sort((a,b)=>a-b).filter(x=>T(x)>=tL-p&&T(x)<=tH+p);
  if(inRange.length<2){
    const rawLo=Math.min(dMin,dMax);const rawHi=Math.max(dMin,dMax);const rawSpan=Math.max(rawHi-rawLo,1);
    const rough=rawSpan/4;const mag=10**Math.floor(Math.log10(rough));
    const step=[1,2,5,10].find(m=>rough<=m*mag)*mag;
    const ticks=[];
    for(let v=Math.ceil(rawLo/step)*step;v<=rawHi+step*0.5;v+=step)ticks.push(Number(v.toFixed(6)));
    if(rawLo<0&&rawHi>0&&!ticks.includes(0))ticks.push(0);
    inRange=ticks.sort((a,b)=>a-b);
  }
  const minGap=span*0.18;
  const thinned=[];
  const has0=inRange.includes(0);
  for(let i=0;i<inRange.length;i++){
    const tv=T(inRange[i]);
    if(thinned.length===0){thinned.push(inRange[i]);continue;}
    const prevTv=T(thinned[thinned.length-1]);
    if(tv-prevTv>=minGap){
      thinned.push(inRange[i]);
    } else if(inRange[i]===0){
      thinned[thinned.length-1]=0;
    }
  }
  if(has0&&!thinned.includes(0)){
    thinned.push(0);
    thinned.sort((a,b)=>a-b);
    const final=[thinned[0]];
    for(let i=1;i<thinned.length;i++){
      const tv=T(thinned[i]);const prevTv=T(final[final.length-1]);
      if(tv-prevTv>=minGap*0.7||thinned[i]===0)final.push(thinned[i]);
    }
    return final;
  }
  return thinned;
}

function analyzeValues(values){
  const sorted=[...values].sort((a,b)=>a-b);
  const n=sorted.length;
  const p005=sorted[Math.floor(n*0.003)];
  const p995=sorted[Math.min(n-1,Math.ceil(n*0.997))];
  const tLo=T(p005);const tHi=T(p995);
  const pad=(tHi-tLo)*0.05;
  return{lo:tLo-pad,hi:tHi+pad,dMin:p005,dMax:p995};
}
function analyzePooledValues(samples,channel){
  const pooled=[];
  for(const s of samples){const vals=s.columns[channel];if(vals)for(let i=0;i<vals.length;i++)pooled.push(vals[i]);}
  if(!pooled.length)return{lo:-1,hi:1,dMin:-1000,dMax:10000};
  pooled.sort((a,b)=>a-b);
  const p005=pooled[Math.floor(pooled.length*0.005)];
  const p995=pooled[Math.ceil(pooled.length*0.995)-1];
  const tLo=T(p005);const tHi=T(p995);
  const pad=(tHi-tLo)*0.06;
  return{lo:tLo-pad,hi:tHi+pad,dMin:p005,dMax:p995};
}
function getRange(values){
  const sorted=[...values].sort((a,b)=>a-b);
  const n=sorted.length;
  const p005=sorted[Math.floor(n*0.003)];
  const p995=sorted[Math.min(n-1,Math.ceil(n*0.997))];
  const tLo=T(p005),tHi=T(p995);
  const pad=(tHi-tLo)*0.05;
  return{lo:tLo-pad,hi:tHi+pad,dMin:p005,dMax:p995};
}
// wider default for scatter/dot plots so high events aren't clipped and there's less dead space
function scatterRange(values){
  const s=[...values].sort((a,b)=>a-b);const n=s.length;
  if(n<3)return{lo:-1,hi:1,dMin:-100,dMax:10000};
  const lo=s[Math.floor(n*0.01)];const hi=s[Math.min(n-1,Math.floor(n*0.9995))];
  const tLo=T(lo),tHi=T(hi);const pad=(tHi-tLo)*0.03||1;
  return{lo:tLo-pad,hi:tHi+pad,dMin:lo,dMax:hi};
}
// resolve manual {min,max} text over an auto range (raw channel units)
function manualBiexp(r,auto){
  if(!auto)return auto;const mn=Number(r.min),mx=Number(r.max);
  if(!r.min.trim()||!r.max.trim()||!isFinite(mn)||!isFinite(mx)||mx<=mn)return auto;
  return{lo:T(mn),hi:T(mx),dMin:mn,dMax:mx};
}
function smoothCounts(counts){
  const kernel=[1,4,6,4,1];
  const radius=2;
  return counts.map((_,i)=>{
    let sum=0;let weight=0;
    for(let k=-radius;k<=radius;k++){
      const idx=i+k;
      if(idx<0||idx>=counts.length)continue;
      const w=kernel[k+radius];
      sum+=counts[idx]*w;
      weight+=w;
    }
    return weight>0?sum/weight:0;
  });
}
function smoothLinePath(points){
  if(points.length===0)return"";
  if(points.length===1)return"M "+points[0].x+","+points[0].y;
  let d="M "+points[0].x+","+points[0].y;
  for(let i=0;i<points.length-1;i++){
    const p1=points[i];
    const p2=points[i+1];
    const mx=(p1.x+p2.x)/2;
    const my=(p1.y+p2.y)/2;
    d+=" Q "+p1.x+","+p1.y+" "+mx+","+my;
  }
  const last=points[points.length-1];
  d+=" Q "+last.x+","+last.y+" "+last.x+","+last.y;
  return d;
}
function smoothAreaPath(points,baseLine){
  if(points.length===0)return"";
  const line=smoothLinePath(points);
  const first=points[0];
  const last=points[points.length-1];
  return line+" L "+last.x+","+baseLine+" L "+first.x+","+baseLine+" Z";
}

// ─── SVG Export helpers ──────────────────────────────
function svgToBlob(svgEl){
  const clone=svgEl.cloneNode(true);
  clone.setAttribute("xmlns","http://www.w3.org/2000/svg");
  clone.style.fontFamily="'IBM Plex Sans', system-ui, sans-serif";
  clone.querySelectorAll("text").forEach(t=>{t.style.fontFamily="'IBM Plex Sans', system-ui, sans-serif";});
  const str=new XMLSerializer().serializeToString(clone);
  return new Blob([str],{type:"image/svg+xml;charset=utf-8"});
}
function downloadBlob(blob,name){
  try{
    const url=URL.createObjectURL(blob);const a=document.createElement("a");
    a.href=url;a.download=name;document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);
  }catch(e){
    try{const url=URL.createObjectURL(blob);window.open(url);}catch(e2){}
  }
}
function downloadDataURL(dataURL,name){
  try{
    const a=document.createElement("a");
    a.href=dataURL;a.download=name;
    document.body.appendChild(a);a.click();document.body.removeChild(a);
  }catch(e){try{window.open(dataURL);}catch(e2){}}
}
function exportSVG(svgEl,name){try{downloadBlob(svgToBlob(svgEl),name+".svg");}catch(e){}}
function exportPNG(svgEl,name,scale=3){
  try{
    const blob=svgToBlob(svgEl);const url=URL.createObjectURL(blob);
    const img=new Image();
    img.onload=()=>{
      try{
        const vb=svgEl.viewBox.baseVal;
        const w=vb.width||svgEl.clientWidth;const h=vb.height||svgEl.clientHeight;
        const canvas=document.createElement("canvas");canvas.width=w*scale;canvas.height=h*scale;
        const ctx=canvas.getContext("2d");ctx.fillStyle="white";ctx.fillRect(0,0,canvas.width,canvas.height);
        ctx.drawImage(img,0,0,canvas.width,canvas.height);
        canvas.toBlob(b=>{if(b)downloadBlob(b,name+".png");},"image/png");
      }catch(e){}
      URL.revokeObjectURL(url);
    };
    img.onerror=()=>{URL.revokeObjectURL(url);};
    img.src=url;
  }catch(e){}
}

// ─── Editable Text ───────────────────────────────────
function EditableText({value,onChange,style,inputStyle}){
  const[editing,setEditing]=useState(false);const[draft,setDraft]=useState(value);const inputRef=useRef(null);
  useEffect(()=>{setDraft(value);},[value]);
  useEffect(()=>{if(editing&&inputRef.current)inputRef.current.select();},[editing]);
  if(editing)return <input ref={inputRef} value={draft} onChange={e=>setDraft(e.target.value)}
    onBlur={()=>{onChange(draft);setEditing(false);}} onKeyDown={e=>{if(e.key==="Enter"){onChange(draft);setEditing(false);}if(e.key==="Escape"){setDraft(value);setEditing(false);}}}
    style={{border:"none",borderBottom:"1.5px solid #3B82F6",outline:"none",background:"transparent",fontFamily:"var(--ff)",padding:"1px 2px",...style,...inputStyle}}/>;
  return <span onClick={()=>setEditing(true)} style={{cursor:"text",borderBottom:"1px dashed transparent",transition:"border-color 0.15s",...style}}
    onMouseEnter={e=>{e.currentTarget.style.borderBottomColor="#CBD5E1";}} onMouseLeave={e=>{e.currentTarget.style.borderBottomColor="transparent";}} title="Click to edit">{value}</span>;
}

// ─── Color Picker ────────────────────────────────────
function normHex(s){let v=(s||"").trim().replace(/^#/,"");if(/^[0-9a-fA-F]{3}$/.test(v))v=v.split("").map(c=>c+c).join("");if(/^[0-9a-fA-F]{6}$/.test(v))return "#"+v.toLowerCase();return null;}
// Strict form used while typing: only a complete 6-digit hex commits, so half-typed input
// (e.g. "#773") is never auto-expanded out from under the caret.
function fullHex(s){const v=(s||"").trim().replace(/^#/,"");return /^[0-9a-fA-F]{6}$/.test(v)?"#"+v.toLowerCase():null;}
function ColorPicker({color,onChange,options=COLOR_OPTIONS,align="left"}){
  const[open,setOpen]=useState(false);const ref=useRef(null);
  const[hex,setHex]=useState(color);
  const editing=useRef(false);
  // Never overwrite the field while the user is typing in it.
  useEffect(()=>{if(!editing.current)setHex(color);},[color]);
  useEffect(()=>{if(!open)return;const c=e=>{if(ref.current&&!ref.current.contains(e.target))setOpen(false);};document.addEventListener("mousedown",c);return()=>document.removeEventListener("mousedown",c);},[open]);
  const safeColor=/^#[0-9a-fA-F]{6}$/.test(color)?color:"#000000";
  const onHexInput=e=>{const raw=e.target.value;setHex(raw);const v=fullHex(raw);if(v)onChange(v);};
  // On blur/Enter, accept 3-digit shorthand; revert to the current color if unparseable.
  const commitHex=()=>{editing.current=false;const v=normHex(hex);if(v){setHex(v);onChange(v);}else setHex(color);};
  return <div ref={ref} style={{position:"relative",display:"inline-block"}}>
    <div onClick={()=>setOpen(!open)} style={{width:18,height:18,borderRadius:5,background:color,cursor:"pointer",border:"2px solid white",boxShadow:"0 0 0 1px #D1D5DB",flexShrink:0}} title="Change color"/>
    {open&&<div style={{position:"absolute",top:24,[align]:-4,zIndex:50,background:"white",borderRadius:10,border:"1px solid #E5E7EB",padding:8,boxShadow:"0 8px 24px rgba(0,0,0,0.12)",width:200}}>
      <div style={{display:"flex",alignItems:"center",gap:6,marginBottom:8}}>
        <input type="color" value={safeColor} onChange={e=>onChange(e.target.value)} title="Color wheel"
          style={{width:34,height:30,padding:0,border:"1px solid #E5E7EB",borderRadius:6,background:"none",cursor:"pointer",flexShrink:0}}/>
        <input value={hex} onChange={onHexInput} maxLength={7} spellCheck={false} placeholder="#RRGGBB"
          onFocus={e=>{editing.current=true;e.target.select();}} onBlur={commitHex}
          onKeyDown={e=>{if(e.key==="Enter"){commitHex();setOpen(false);}else if(e.key==="Escape"){editing.current=false;setHex(color);}}}
          style={{flex:1,minWidth:0,padding:"6px 8px",borderRadius:6,border:"1px solid "+(normHex(hex)?"#D1D5DB":"#FCA5A5"),fontSize:12.5,fontFamily:"ui-monospace,Menlo,monospace",color:"#111827",boxSizing:"border-box"}}/>
      </div>
      <div style={{display:"grid",gridTemplateColumns:"repeat(7,1fr)",gap:5,maxHeight:170,overflowY:"auto"}}>
        {options.map(c=><div key={c} onClick={()=>{onChange(c);setOpen(false);}} style={{width:22,height:22,borderRadius:5,background:c,cursor:"pointer",border:c.toLowerCase()===(color||"").toLowerCase()?"2.5px solid #111":"2px solid transparent",transition:"transform 0.1s"}}
          onMouseEnter={e=>{e.currentTarget.style.transform="scale(1.2)";}} onMouseLeave={e=>{e.currentTarget.style.transform="scale(1)";}}/>)}
      </div>
    </div>}
  </div>;
}

// ─── Single Histogram ────────────────────────────────
function Histogram({values,name,color,xLabel,yLabel,gateValue,onGateChange,onNameChange,xDomain,svgRef:externalRef,gateLabel,gateColor,yMode="count",showGate=true,showPct=true}){
  const internalRef=useRef(null);
  const svgRef=externalRef||internalRef;
  const dragRef=useRef(false);

  const{lo,hi,dMin,dMax}=useMemo(()=>xDomain||analyzeValues(values),[xDomain,values]);
  const{counts,centers}=useMemo(()=>computeHist(values,N_BINS,lo,hi),[values,lo,hi]);
  const displayCounts=useMemo(()=>smoothCounts(counts),[counts]);
  const maxC=useMemo(()=>{let m=1;for(let i=0;i<displayCounts.length;i++)if(displayCounts[i]>m)m=displayCounts[i];return m;},[displayCounts]);

  const xS=useCallback(tv=>PL.ml+((tv-lo)/(hi-lo))*PW,[lo,hi]);
  const yS=useCallback(c=>PL.mt+PH-(c/maxC)*PH,[maxC]);

  const tGate=T(gateValue);
  const gateX=Math.max(PL.ml,Math.min(PL.ml+PW,xS(tGate)));
  const pePct=useMemo(()=>{let c=0;for(let i=0;i<values.length;i++)if(values[i]>=gateValue)c++;return((c/values.length)*100).toFixed(1);},[values,gateValue]);
  const gmfiAll=useMemo(()=>{let logSum=0,logN=0;for(let i=0;i<values.length;i++)if(values[i]>0){logSum+=Math.log(values[i]);logN++;}return logN>0?Math.round(Math.exp(logSum/logN)):0;},[values]);

  const{mainPath,posPath}=useMemo(()=>{
    const baseLine=PL.mt+PH;
    const points=centers.map((center,i)=>({x:xS(center),y:yS(displayCounts[i])}));
    const mp=smoothAreaPath(points,baseLine);
    let pp="";const tg=T(gateValue);const gi=centers.findIndex(c=>c>=tg);
    if(gi>=0&&gi<centers.length){
      const gx=Math.max(PL.ml,Math.min(PL.ml+PW,xS(tg)));
      let sY=baseLine;
      if(gi>0){
        const f=(tg-centers[gi-1])/(centers[gi]-centers[gi-1]);
        sY=yS(displayCounts[gi-1]+f*(displayCounts[gi]-displayCounts[gi-1]));
      }else sY=yS(displayCounts[0]);
      const tail=[{x:gx,y:sY},...points.slice(gi).filter(pt=>pt.x>gx)];
      pp=smoothAreaPath(tail,baseLine);
    }
    return{mainPath:mp,posPath:pp};
  },[centers,displayCounts,gateValue,xS,yS]);

  const xticks=useMemo(()=>axisTicks(dMin,dMax,PW),[dMin,dMax]);
  const yTicks=useMemo(()=>{
    if(yMode==="pct")return[0,0.25,0.5,0.75,1].map(f=>({y:PL.mt+PH-f*PH,label:Math.round(f*100)}));
    const n=4;const s=Math.ceil(maxC/n);
    return Array.from({length:n+1},(_,i)=>i*s).filter(v=>v<=maxC*1.05).map(v=>({y:yS(v),label:v.toLocaleString()}));
  },[maxC,yMode,yS]);
  const effYLabel=yMode==="pct"?"% of max":yLabel;

  const toVal=useCallback(e=>{if(!svgRef.current)return gateValue;const r=svgRef.current.getBoundingClientRect();const sx=((e.clientX-r.left)/r.width)*PL.w;const tv=lo+((sx-PL.ml)/PW)*(hi-lo);return Math.round(invT(Math.max(lo,Math.min(hi,tv))));},[lo,hi,gateValue,svgRef]);
  const onDown=e=>{dragRef.current=true;e.preventDefault();};
  const onMove=useCallback(e=>{if(dragRef.current)onGateChange(toVal(e));},[toVal,onGateChange]);
  const onUp=()=>{dragRef.current=false;};

  const flipAnn=gateX+72>PL.ml+PW;

  return(
    <div style={{background:"white",borderRadius:10,border:"1px solid #e2e5ea",padding:"10px 6px 6px",boxShadow:"0 1px 4px rgba(0,0,0,0.04)"}}>
      <div style={{display:"flex",alignItems:"center",justifyContent:"center",gap:6,marginBottom:2,minHeight:20}}>
        <EditableText value={name} onChange={onNameChange} style={{fontSize:12.5,fontWeight:500,color:"#111827"}} inputStyle={{fontSize:12.5,fontWeight:500,textAlign:"center",width:"70%"}}/>
        <button onClick={()=>svgRef.current&&exportSVG(svgRef.current,name.replace(/\s+/g,"_"))} title="Export SVG"
          style={{background:"none",border:"none",cursor:"pointer",padding:2,color:"#CBD5E1",display:"flex"}}
          onMouseEnter={e=>{e.currentTarget.style.color="#6B7280";}} onMouseLeave={e=>{e.currentTarget.style.color="#CBD5E1";}}>
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
        </button>
      </div>
      <div style={{textAlign:"center",fontSize:9.5,color:"#6B7280",fontWeight:600,marginTop:-1,marginBottom:2}}>{"gMFI "+gmfiAll.toLocaleString()}</div>
      <svg ref={svgRef} viewBox={"0 0 "+PL.w+" "+PL.h} style={{width:"100%",display:"block",userSelect:"none"}} onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={onUp}>
        {yTicks.map((yt,i)=><line key={i} x1={PL.ml} x2={PL.ml+PW} y1={yt.y} y2={yt.y} stroke="#f0f1f3" strokeWidth={0.8}/>)}
        <path d={mainPath} fill={color+"28"} stroke={color} strokeWidth={1} strokeLinejoin="round"/>
        {showGate&&posPath&&<path d={posPath} fill={gateColor+"22"} stroke="none"/>}
        {showGate&&<>
        <line x1={gateX} x2={gateX} y1={PL.mt} y2={PL.mt+PH} stroke={gateColor} strokeWidth={0.9} strokeDasharray="4,3"/>
        <rect x={gateX-12} y={PL.mt-4} width={24} height={PH+8} fill="transparent" style={{cursor:"ew-resize"}} onMouseDown={onDown}/>
        </>}
        {showPct&&<g transform={"translate("+(flipAnn?gateX-68:gateX+10)+","+(PL.mt+14)+")"}>
          <rect x={-5} y={-13} width={62} height={36} rx={5} fill="white" fillOpacity={0.92} stroke={gateColor+"55"} strokeWidth={0.8}/>
          <text fontSize="8.5" fill="#9CA3AF" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{gateLabel}</text>
          <text y={16} fontSize="13" fontWeight="800" fill={gateColor} style={{fontFamily:"var(--ff)"}}>{pePct}%</text>
        </g>}
        <line x1={PL.ml} x2={PL.ml+PW} y1={PL.mt+PH} y2={PL.mt+PH} stroke="#111111" strokeWidth={1.2}/>
        {xticks.map(({v,label})=>{const x=xS(T(v));return <g key={v}><line x1={x} x2={x} y1={PL.mt+PH} y2={PL.mt+PH+(label?6:3)} stroke={label?"#111111":"#9CA3AF"}/>{label&&<text x={x} y={PL.mt+PH+17} textAnchor="middle" fontSize="9.5" fill="#4B5563" style={{fontFamily:"var(--ff)"}}>{fmtTick(v)}</text>}</g>;})}
        <text x={PL.ml+PW/2} y={PL.h-4} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{xLabel}</text>
        <line x1={PL.ml} x2={PL.ml} y1={PL.mt} y2={PL.mt+PH} stroke="#4B5563" strokeWidth={1}/>
        {yTicks.map((yt,i)=><g key={i}><line x1={PL.ml-4} x2={PL.ml} y1={yt.y} y2={yt.y} stroke="#4B5563"/><text x={PL.ml-7} y={yt.y+3} textAnchor="end" fontSize="8" fill="#6B7280" style={{fontFamily:"var(--ff)"}}>{yt.label}</text></g>)}
        <text transform={"translate(13,"+(PL.mt+PH/2)+") rotate(-90)"} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{effYLabel}</text>
        <text x={PL.ml+PW-2} y={PL.mt+PH-6} textAnchor="end" fontSize="8" fill="#9CA3AF" style={{fontFamily:"var(--ff)"}}>{"n = "+values.length.toLocaleString()}</text>
        <text x={PL.ml+PW-2} y={PL.mt+PH-18} textAnchor="end" fontSize="8" fill="#9CA3AF" style={{fontFamily:"var(--ff)"}}>{"gMFI = "+gmfiAll.toLocaleString()}</text>
      </svg>
    </div>
  );
}

// ─── Overlay Histogram (NEW) ─────────────────────────
// Superimposes several samples' density curves on one shared axis.
function OverlayHistogram({samples,colors,channel,xLabel,yLabel,gateValue,onGateChange,xDomain,gateLabel,normalize,gateColor,showGate=true,showPct=true,onToggleGate,onTogglePct,pctPos="left",setPctPos}){
  const svgRef=useRef(null);
  const dragRef=useRef(false);

  const{lo,hi,dMin,dMax}=useMemo(()=>xDomain||analyzePooledValues(samples,channel),[xDomain,samples,channel]);

  const series=useMemo(()=>samples.map((s,i)=>{
    const vals=s.columns[channel]||[];
    const{counts,centers}=computeHist(vals,N_BINS,lo,hi);
    const displayCounts=smoothCounts(counts);
    let mx=1;for(let k=0;k<displayCounts.length;k++)if(displayCounts[k]>mx)mx=displayCounts[k];
    let pos=0;for(let k=0;k<vals.length;k++)if(vals[k]>=gateValue)pos++;
    let logSum=0,logN=0;for(let k=0;k<vals.length;k++)if(vals[k]>0){logSum+=Math.log(vals[k]);logN++;}
    const gmfi=logN>0?Math.round(Math.exp(logSum/logN)):0;
    return{name:s.name,color:colors[i]||PALETTE[i%PALETTE.length],counts:displayCounts,centers,maxC:mx,n:vals.length,pct:vals.length?((pos/vals.length)*100).toFixed(1):"0.0",gmfi};
  }),[samples,colors,channel,lo,hi,gateValue]);

  const globalMax=useMemo(()=>{let m=1;for(const s of series)if(s.maxC>m)m=s.maxC;return m;},[series]);

  const xS=useCallback(tv=>PL.ml+((tv-lo)/(hi-lo))*PW,[lo,hi]);
  const yS=useCallback((c,ref)=>PL.mt+PH-(c/ref)*PH,[]);

  const paths=useMemo(()=>{
    const baseLine=PL.mt+PH;
    return series.map(s=>{
      const ref=normalize?s.maxC:globalMax;
      const points=s.centers.map((center,i)=>({x:xS(center),y:yS(s.counts[i],ref)}));
      return{line:smoothLinePath(points),area:smoothAreaPath(points,baseLine),color:s.color};
    });
  },[series,normalize,globalMax,xS,yS]);

  // right-shoulder anchor (~60% of peak height) for peak-placed % labels, per curve
  const peakLabels=useMemo(()=>series.map(s=>{
    let pk=0;for(let j=1;j<s.counts.length;j++)if(s.counts[j]>s.counts[pk])pk=j;
    const ref=normalize?s.maxC:globalMax;
    const thr=s.counts[pk]*0.6;
    let ri=pk;for(let j=pk+1;j<s.counts.length;j++){ri=j;if(s.counts[j]<=thr)break;}
    const x=Math.min(xS(s.centers[ri])+3,PL.ml+PW-30);
    const y=Math.max(PL.mt+8,yS(s.counts[ri],ref)-3);
    return{x,y,pct:s.pct,color:s.color};
  }),[series,normalize,globalMax,xS,yS]);

  const xticks=useMemo(()=>axisTicks(dMin,dMax,PW),[dMin,dMax]);
  const yTicks=useMemo(()=>{
    if(normalize)return[0,0.25,0.5,0.75,1];
    const n=4;const s=Math.ceil(globalMax/n);return Array.from({length:n+1},(_,i)=>i*s).filter(v=>v<=globalMax*1.05);
  },[normalize,globalMax]);
  const yTickY=useCallback(yt=>normalize?PL.mt+PH-yt*PH:yS(yt,globalMax),[normalize,globalMax,yS]);
  const yTickLabel=useCallback(yt=>normalize?Math.round(yt*100):yt,[normalize]);

  const tGate=T(gateValue);
  const gateX=Math.max(PL.ml,Math.min(PL.ml+PW,xS(tGate)));

  const toVal=useCallback(e=>{if(!svgRef.current)return gateValue;const r=svgRef.current.getBoundingClientRect();const sx=((e.clientX-r.left)/r.width)*PL.w;const tv=lo+((sx-PL.ml)/PW)*(hi-lo);return Math.round(invT(Math.max(lo,Math.min(hi,tv))));},[lo,hi,gateValue]);
  const onDown=e=>{dragRef.current=true;e.preventDefault();};
  const onMove=useCallback(e=>{if(dragRef.current)onGateChange(toVal(e));},[toVal,onGateChange]);
  const onUp=()=>{dragRef.current=false;};

  const effYLabel=normalize?"% of max":yLabel;
  const OW=PL.w+170; // extra width on the right for the legend (swatch + name)
  const flipAnn=gateX+72>PL.ml+PW;

  return(
    <div style={{background:"white",borderRadius:12,border:"1px solid #e2e5ea",padding:"14px 12px",boxShadow:"0 1px 4px rgba(0,0,0,0.04)"}}>
      <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",marginBottom:6,padding:"0 4px"}}>
        <span style={{fontSize:13,fontWeight:700,color:"#111827"}}>Overlay — {channel} · {series.length} samples</span>
        <div style={{display:"flex",gap:6,flexWrap:"wrap",justifyContent:"flex-end"}}>
          {onToggleGate&&<button onClick={onToggleGate} title="Show/hide the gate line" style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(showGate?"#BFDBFE":"#E5E7EB"),background:showGate?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:showGate?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>Gate</button>}
          {onTogglePct&&<button onClick={onTogglePct} title="Show/hide % positive labels" style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(showPct?"#BFDBFE":"#E5E7EB"),background:showPct?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:showPct?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>% positive</button>}
          {showPct&&<div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <button onClick={()=>setPctPos("left")} title="% labels in the legend" style={{padding:"4px 9px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:pctPos==="left"?"#EFF6FF":"white",color:pctPos==="left"?"#3B82F6":"#9CA3AF"}}>Legend</button>
            <button onClick={()=>setPctPos("peak")} title="% labels on the curves, right of each peak" style={{padding:"4px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:pctPos==="peak"?"#EFF6FF":"white",color:pctPos==="peak"?"#3B82F6":"#9CA3AF"}}>Peak</button>
          </div>}
          <button onClick={()=>svgRef.current&&exportSVG(svgRef.current,"overlay")} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>SVG</button>
          <button onClick={()=>svgRef.current&&exportPNG(svgRef.current,"overlay")} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>PNG</button>
        </div>
      </div>
      <svg ref={svgRef} viewBox={"0 0 "+OW+" "+PL.h} style={{width:"100%",display:"block",userSelect:"none"}} onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={onUp}>
        {yTicks.map((yt,i)=><line key={i} x1={PL.ml} x2={PL.ml+PW} y1={yTickY(yt)} y2={yTickY(yt)} stroke="#f0f1f3" strokeWidth={0.8}/>)}
        {/* fills first (behind), then strokes on top for clean overlapping lines */}
        {paths.map((p,i)=><path key={"a"+i} d={p.area} fill={p.color+"1E"} stroke="none"/>)}
        {paths.map((p,i)=><path key={"l"+i} d={p.line} fill="none" stroke={p.color} strokeWidth={1.1} strokeLinejoin="round"/>)}
        {showGate&&<>
        <line x1={gateX} x2={gateX} y1={PL.mt} y2={PL.mt+PH} stroke={gateColor} strokeWidth={0.9} strokeDasharray="4,3"/>
        <rect x={gateX-12} y={PL.mt-4} width={24} height={PH+8} fill="transparent" style={{cursor:"ew-resize"}} onMouseDown={onDown}/>
        <text x={flipAnn?gateX-8:gateX+8} y={PL.mt+9} textAnchor={flipAnn?"end":"start"} fontSize="8.5" fill="#9CA3AF" fontWeight="600" style={{fontFamily:"var(--ff)"}}>{gateLabel}</text>
        </>}
        {/* axes */}
        <line x1={PL.ml} x2={PL.ml+PW} y1={PL.mt+PH} y2={PL.mt+PH} stroke="#111111" strokeWidth={1.2}/>
        {xticks.map(({v,label})=>{const x=xS(T(v));return <g key={v}><line x1={x} x2={x} y1={PL.mt+PH} y2={PL.mt+PH+(label?6:3)} stroke={label?"#111111":"#9CA3AF"}/>{label&&<text x={x} y={PL.mt+PH+17} textAnchor="middle" fontSize="9.5" fill="#4B5563" style={{fontFamily:"var(--ff)"}}>{fmtTick(v)}</text>}</g>;})}
        <text x={PL.ml+PW/2} y={PL.h-4} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{xLabel}</text>
        <line x1={PL.ml} x2={PL.ml} y1={PL.mt} y2={PL.mt+PH} stroke="#4B5563" strokeWidth={1}/>
        {yTicks.map((yt,i)=><g key={i}><line x1={PL.ml-4} x2={PL.ml} y1={yTickY(yt)} y2={yTickY(yt)} stroke="#4B5563"/><text x={PL.ml-7} y={yTickY(yt)+3} textAnchor="end" fontSize="8" fill="#6B7280" style={{fontFamily:"var(--ff)"}}>{yTickLabel(yt)}</text></g>)}
        <text transform={"translate(13,"+(PL.mt+PH/2)+") rotate(-90)"} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{effYLabel}</text>
        {/* peak-placed % labels on the curves */}
        {showPct&&pctPos==="peak"&&peakLabels.map((p,i)=><text key={"pk"+i} x={p.x} y={p.y} textAnchor="start" fontSize="8.5" fontWeight="700" fill={p.color} style={{fontFamily:"var(--ff)"}}>{p.pct}%</text>)}
        {/* legend on the right — each name grouped with its %, blocks evenly spaced */}
        {(()=>{
          const bh=Math.min(28,(PH-6)/series.length);
          const startY=PL.mt+Math.max(4,(PH-series.length*bh)/2);
          const pctDy=Math.min(11,bh*0.42);
          return series.map((s,i)=>(
            <g key={"leg"+i} transform={"translate("+(PL.ml+PW+16)+","+(startY+i*bh)+")"}>
              {/* color key — always shown so overlapping curves are identifiable */}
              <rect x={0} y={-8} width={9} height={9} rx={2} fill={s.color}/>
              <text x={15} y={0} fontSize="9" fill="#374151" fontWeight="600" style={{fontFamily:"var(--ff)"}}>{s.name.length>20?s.name.slice(0,19)+"…":s.name}</text>
              {showPct&&pctPos==="left"&&<text x={15} y={pctDy} fontSize="8" fill={gateColor} fontWeight="500" style={{fontFamily:"var(--ff)"}}>{s.pct+"% "+gateLabel}</text>}
            </g>
          ));
        })()}
      </svg>
    </div>
  );
}

// ─── Ridge Plot ──────────────────────────────────────
function RidgePlot({samples,colors,channel,gateValue,onGateChange,xLabel,xDomain,gateLabel,gateColor,showGate=true,onToggleGate,showPct=true,onTogglePct,pctPos="left",setPctPos}){
  const ridgeRef=useRef(null);
  const dragRef=useRef(false);
  const[labelMode,setLabelMode]=useState("side"); // "side" (left gutter) | "inline" (compact, on plot) | "legend" (key on right)
  const[fontScale,setFontScale]=useState(1);
  const[overlap,setOverlap]=useState(0.55); // 0 = staggered/separated lanes … 0.7 = tight ridge
  const[legendStats,setLegendStats]=useState(true); // show gMFI/% inside the legend key
  const[statsOnPlot,setStatsOnPlot]=useState(false); // draw gMFI/% as a white badge on each row instead

  const ROW_H=92;
  // 0 = fully separated lanes (staggered, nothing can occlude), up to 0.7 = tightly stacked ridge
  const STEP=ROW_H*(1-overlap);
  // In legend mode, reserve a top band the height of the key so the curves shift down and the
  // key never overlaps a plot. The band grows with the number of samples and the font size.
  const legLineH=Math.max(8,12.5*fontScale);
  // wrap the key into columns (max 4 rows each) so many samples stay compact instead of one tall stack
  const legCols=labelMode==="legend"?Math.max(1,Math.ceil(samples.length/4)):1;
  const legPerCol=Math.max(1,Math.ceil(samples.length/legCols));
  const legTopPad=labelMode==="legend"?legPerCol*legLineH+16:0;
  const MT=20+legTopPad,MB=64;
  const PLOT_L=10,PLOT_R=10;
  const PLOT_W=500;
  const n=samples.length;
  const svgH=MT+ROW_H+STEP*(n-1)+MB;
  const lastRowBase=MT+ROW_H+STEP*(n-1);
  const maxNameLen=Math.max(1,...samples.map(s=>s.name.length));
  const LABEL_W=(labelMode==="inline"||labelMode==="legend")?8:Math.max(maxNameLen*7.5*fontScale,66)+22;
  const vbX=-LABEL_W;
  const vbW=LABEL_W+PLOT_L+PLOT_W+PLOT_R; // legend overlaps the top-right corner, no extra width

  const sharedRange=useMemo(()=>xDomain||analyzePooledValues(samples,channel),[xDomain,samples,channel]);
  const xS=useCallback(tv=>PLOT_L+((tv-sharedRange.lo)/(sharedRange.hi-sharedRange.lo))*PLOT_W,[sharedRange]);
  const xticks=useMemo(()=>axisTicks(sharedRange.dMin,sharedRange.dMax,PLOT_W),[sharedRange]);

  const histData=useMemo(()=>{
    return samples.map(s=>{
      const vals=s.columns[channel]||[];
      const{counts,centers}=computeHist(vals,N_BINS,sharedRange.lo,sharedRange.hi);
      const displayCounts=smoothCounts(counts);
      let maxC=1;for(let i=0;i<displayCounts.length;i++)if(displayCounts[i]>maxC)maxC=displayCounts[i];
      let posCount=0;for(let i=0;i<vals.length;i++)if(vals[i]>=gateValue)posCount++;
      const pct=vals.length>0?((posCount/vals.length)*100).toFixed(1):"0.0";
      let logSum=0,logN=0;for(let i=0;i<vals.length;i++)if(vals[i]>0){logSum+=Math.log(vals[i]);logN++;}
      const gmfi=logN>0?Math.round(Math.exp(logSum/logN)):0;
      return{counts:displayCounts,centers,maxC,pct,n:vals.length,gmfi};
    });
  },[samples,channel,sharedRange,gateValue]);

  const tGate=T(gateValue);
  const gateX=Math.max(PLOT_L,Math.min(PLOT_L+PLOT_W,xS(tGate)));

  const toVal=useCallback(e=>{
    if(!ridgeRef.current)return gateValue;const r=ridgeRef.current.getBoundingClientRect();
    const pxPerUnit=r.width/vbW;
    const svgX=((e.clientX-r.left)/pxPerUnit)+vbX;
    const tv=sharedRange.lo+((svgX-PLOT_L)/PLOT_W)*(sharedRange.hi-sharedRange.lo);
    return Math.round(invT(Math.max(sharedRange.lo,Math.min(sharedRange.hi,tv))));
  },[sharedRange,gateValue,vbW,vbX]);
  const onDown=e=>{dragRef.current=true;e.preventDefault();};
  const onMove=useCallback(e=>{if(dragRef.current)onGateChange(toVal(e));},[toVal,onGateChange]);
  const onUp=()=>{dragRef.current=false;};

  const rows=samples.map((s,i)=>{
    const hd=histData[i];
    const yOff=MT+i*STEP;
    const localYS=c=>yOff+ROW_H-(c/hd.maxC)*ROW_H*0.85;
    const color=colors[i]||PALETTE[i%PALETTE.length];
    const baseLine=yOff+ROW_H;
    const points=hd.centers.map((center,j)=>({x:xS(center),y:localYS(hd.counts[j])}));
    const path=smoothAreaPath(points,baseLine);
    const line=smoothLinePath(points);
    let pk=0;for(let j=1;j<hd.counts.length;j++)if(hd.counts[j]>hd.counts[pk])pk=j;
    // anchor the peak-placed % label on the curve's right shoulder (~60% of peak height)
    const thr=hd.counts[pk]*0.6;
    let ri=pk;for(let j=pk+1;j<hd.counts.length;j++){ri=j;if(hd.counts[j]<=thr)break;}
    const labelX=Math.min(xS(hd.centers[ri])+3,PLOT_L+PLOT_W-30);
    const labelY=Math.min(baseLine-4,localYS(hd.counts[ri])-3);
    return{yOff,baseLine,path,line,color,name:s.name,pct:hd.pct,gmfi:hd.gmfi,labelX,labelY};
  });

  return (
    <div style={{background:"white",borderRadius:12,border:"1px solid #e2e5ea",padding:"16px 12px",boxShadow:"0 1px 4px rgba(0,0,0,0.04)"}}>
      <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",marginBottom:8,padding:"0 4px"}}>
        <span style={{fontSize:13,fontWeight:700,color:"#111827"}}>Ridge Plot — {channel}</span>
        <div style={{display:"flex",gap:6}}>
          {onToggleGate&&<button onClick={onToggleGate} title="Show/hide the gate line" style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(showGate?"#BFDBFE":"#E5E7EB"),background:showGate?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:showGate?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>Gate</button>}
          {onTogglePct&&<button onClick={onTogglePct} title="Show/hide % positive per row" style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(showPct?"#BFDBFE":"#E5E7EB"),background:showPct?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:showPct?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>% positive</button>}
          {showPct&&<div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <button onClick={()=>setPctPos("left")} title="% label by the sample name" style={{padding:"4px 9px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:pctPos==="left"?"#EFF6FF":"white",color:pctPos==="left"?"#3B82F6":"#9CA3AF"}}>Left</button>
            <button onClick={()=>setPctPos("peak")} title="% label on the histogram, right of the peak" style={{padding:"4px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:pctPos==="peak"?"#EFF6FF":"white",color:pctPos==="peak"?"#3B82F6":"#9CA3AF"}}>Peak</button>
          </div>}
          <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <button onClick={()=>setLabelMode("side")} title="Labels in a left gutter" style={{padding:"4px 9px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:labelMode==="side"?"#EFF6FF":"white",color:labelMode==="side"?"#3B82F6":"#9CA3AF"}}>Side</button>
            <button onClick={()=>setLabelMode("inline")} title="Compact: labels above each baseline (fits slides)" style={{padding:"4px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:labelMode==="inline"?"#EFF6FF":"white",color:labelMode==="inline"?"#3B82F6":"#9CA3AF"}}>Compact</button>
            <button onClick={()=>setLabelMode("legend")} title="Separate key on the right, nothing on the plot" style={{padding:"4px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:labelMode==="legend"?"#EFF6FF":"white",color:labelMode==="legend"?"#3B82F6":"#9CA3AF"}}>Legend</button>
          </div>
          {labelMode==="legend"&&<button onClick={()=>setLegendStats(v=>!v)} title="Show gMFI / % inside the key, or names only (stats stay in the table below)"
            style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(legendStats?"#BFDBFE":"#E5E7EB"),background:legendStats?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:legendStats?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>Key stats</button>}
          {labelMode==="legend"&&<button onClick={()=>setStatsOnPlot(v=>!v)} title="Draw gMFI / % as a small white badge on each histogram row"
            style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(statsOnPlot?"#BFDBFE":"#E5E7EB"),background:statsOnPlot?"#EFF6FF":"white",fontSize:11,fontWeight:600,color:statsOnPlot?"#3B82F6":"#9CA3AF",cursor:"pointer",fontFamily:"var(--ff)"}}>Stats on plot</button>}
          <div style={{display:"flex",alignItems:"center",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}} title="Row overlap — 0% fully separates the rows so curves can't cover each other">
            <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",padding:"0 6px"}}>Overlap</span>
            <button onClick={()=>setOverlap(o=>Math.max(0,+(o-0.1).toFixed(2)))} style={{padding:"3px 8px",border:"none",borderLeft:"1px solid #E5E7EB",background:"white",fontSize:12,fontWeight:700,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>−</button>
            <span style={{fontSize:10,fontWeight:600,color:"#9CA3AF",minWidth:30,textAlign:"center"}}>{Math.round(overlap*100)}%</span>
            <button onClick={()=>setOverlap(o=>Math.min(0.7,+(o+0.1).toFixed(2)))} style={{padding:"3px 8px",border:"none",borderLeft:"1px solid #E5E7EB",background:"white",fontSize:12,fontWeight:700,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>+</button>
          </div>
          <div style={{display:"flex",alignItems:"center",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}} title="Label font size">
            <button onClick={()=>setFontScale(s=>Math.max(0.35,+(s-0.05).toFixed(2)))} style={{padding:"3px 8px",border:"none",background:"white",fontSize:10,fontWeight:700,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>A−</button>
            <span style={{fontSize:10,fontWeight:600,color:"#9CA3AF",minWidth:26,textAlign:"center"}}>{Math.round(fontScale*100)}%</span>
            <button onClick={()=>setFontScale(s=>Math.min(2,+(s+0.05).toFixed(2)))} style={{padding:"3px 8px",border:"none",borderLeft:"1px solid #E5E7EB",background:"white",fontSize:13,fontWeight:700,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>A+</button>
          </div>
          <button onClick={()=>ridgeRef.current&&exportSVG(ridgeRef.current,"ridge_plot")} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>SVG</button>
          <button onClick={()=>ridgeRef.current&&exportPNG(ridgeRef.current,"ridge_plot")} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>PNG</button>
        </div>
      </div>
      <svg ref={ridgeRef} viewBox={vbX+" 0 "+vbW+" "+svgH} style={{width:"100%",display:"block",userSelect:"none",background:"white",overflow:"visible"}}
        onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={onUp}>
        {rows.map((r,i)=>{
          const bl=r.yOff+ROW_H; // this ridge's own baseline
          return(
          <g key={i}>
            <path d={r.path} fill="white" stroke="none"/>
            <path d={r.path} fill={r.color+"25"} stroke="none"/>
            {/* each row gets its own black x-axis baseline */}
            <line x1={PLOT_L} x2={PLOT_L+PLOT_W} y1={bl} y2={bl} stroke="#111111" strokeWidth={1}/>
            <path d={r.line} fill="none" stroke={r.color} strokeWidth={1.1} strokeLinejoin="round"/>
          </g>
        );})}
        {showGate&&<>
        <line x1={gateX} x2={gateX} y1={MT} y2={MT+ROW_H+STEP*(n-1)} stroke={gateColor} strokeWidth={0.9} strokeDasharray="4,3"/>
        <rect x={gateX-12} y={MT-4} width={24} height={ROW_H+STEP*(n-1)+8} fill="transparent" style={{cursor:"ew-resize"}} onMouseDown={onDown}/>
        </>}
        <line x1={PLOT_L} x2={PLOT_L+PLOT_W} y1={lastRowBase} y2={lastRowBase} stroke="#111111" strokeWidth={1.2}/>
        {xticks.map(({v,label})=>{const x=xS(T(v));return <g key={v}><line x1={x} x2={x} y1={lastRowBase} y2={lastRowBase+(label?6:3)} stroke={label?"#111111":"#9CA3AF"}/>{label&&<text x={x} y={lastRowBase+11} textAnchor="middle" dominantBaseline="hanging" fontSize="10" fill="#4B5563" style={{fontFamily:"var(--ff)"}}>{fmtTick(v)}</text>}</g>;})}
        {/* legend mode: a compact key overlapping the top-right corner, nothing on the plot */}
        {labelMode==="legend"&&(()=>{
          const leftPct=showPct&&pctPos==="left";
          const fs=9.5*fontScale,lineH=legLineH;
          const cw=5.3*fontScale;                       // approx glyph width at this size
          const sw=Math.max(5,7*fontScale);             // swatch size
          const gap=sw+4,padX=6;                        // swatch→text gap, column padding
          const raw=rows.map(r=>({color:r.color,name:r.name,
            sub:legendStats?(" · gMFI "+r.gmfi.toLocaleString()+(leftPct?" · "+r.pct+"%":"")):""}));
          // Truncate names so the whole key fits inside the plot width, never wider than the figure.
          const perColChars=Math.floor((PLOT_W/legCols-gap-padX*2)/cw);
          const entries=raw.map(e=>{
            const total=e.name.length+e.sub.length;
            if(total<=perColChars)return e;
            const keep=Math.max(5,perColChars-e.sub.length);
            return{...e,name:e.name.length>keep?e.name.slice(0,keep-1)+"…":e.name};
          });
          // size each column to its own longest entry, not the global longest
          const colWs=[];
          for(let c=0;c<legCols;c++){
            let m=0;
            for(let r2=0;r2<legPerCol;r2++){const i=c*legPerCol+r2;if(i<entries.length)m=Math.max(m,(entries[i].name+entries[i].sub).length);}
            colWs.push(gap+m*cw+padX*2);
          }
          const colX=[];let acc=0;for(let c=0;c<legCols;c++){colX.push(acc);acc+=colWs[c];}
          const boxW=acc+4,boxH=legPerCol*lineH+7;
          const bx=Math.max(PLOT_L,PLOT_L+PLOT_W-boxW),by=4;
          return(
            <g>
              <rect x={bx} y={by} width={boxW} height={boxH} rx={4} fill="white" fillOpacity={0.85} stroke="#E5E7EB" strokeWidth={0.7}/>
              {entries.map((e,i)=>{
                const col=Math.floor(i/legPerCol),row=i%legPerCol;
                return(
                  <g key={"leg"+i} transform={"translate("+(bx+padX+colX[col])+","+(by+lineH*0.78+row*lineH)+")"}>
                    <rect x={0} y={-sw*0.85} width={sw} height={sw} rx={1.5} fill={e.color}/>
                    <text x={gap} y={0} fontSize={fs} fontWeight="600" fill="#374151" style={{fontFamily:"var(--ff)"}}>{e.name}<tspan fill="#9CA3AF" fontWeight="500">{e.sub}</tspan></text>
                  </g>
                );
              })}
            </g>
          );
        })()}
        {/* stats as a white badge on each row, so the key can stay names-only */}
        {labelMode==="legend"&&statsOnPlot&&rows.map((r,i)=>{
          const fs=Math.max(5,8.5*fontScale),cw=4.9*fontScale;
          const txt="gMFI "+r.gmfi.toLocaleString()+(showPct?"  ·  "+r.pct+"%":"");
          const w=txt.length*cw+9,h=fs+6;
          const x=PLOT_L+PLOT_W-w-3,y=r.yOff+ROW_H-h-3;
          return(
            <g key={"badge"+i} style={{pointerEvents:"none"}}>
              <rect x={x} y={y} width={w} height={h} rx={3} fill="white" fillOpacity={0.92} stroke={r.color} strokeWidth={0.7}/>
              <text x={x+w/2} y={y+h*0.72} textAnchor="middle" fontSize={fs} fontWeight="600" fill="#374151" style={{fontFamily:"var(--ff)"}}>{txt}</text>
            </g>
          );
        })}
        {/* labels drawn last so they sit on top of every curve */}
        {labelMode!=="legend"&&rows.map((r,i)=>{
          const bl=r.yOff+ROW_H;const leftPct=showPct&&pctPos==="left";const lh=10*fontScale;
          if(labelMode==="side")return(
            <g key={"lab"+i}>
              <text x={-8} y={leftPct?bl+1-lh:bl-lh*0.5} textAnchor="end" fontSize={10.5*fontScale} fontWeight="500" fill="#374151" style={{fontFamily:"var(--ff)"}}>{r.name}</text>
              <text x={-8} y={leftPct?bl+1:bl+lh*0.5} textAnchor="end" fontSize={8.5*fontScale} fontWeight="600" fill="#6B7280" style={{fontFamily:"var(--ff)"}}>{"gMFI "+r.gmfi.toLocaleString()}</text>
              {leftPct&&<text x={-8} y={bl+1+lh} textAnchor="end" fontSize={9*fontScale} fontWeight="700" fill={gateColor} style={{fontFamily:"var(--ff)"}}>{r.pct}%</text>}
            </g>
          );
          const gmfiTxt="gMFI "+r.gmfi.toLocaleString()+(leftPct?"  ·  "+r.pct+"%":"");
          return(
            <g key={"lab"+i} style={{pointerEvents:"none"}}>
              <text x={PLOT_L+3} y={bl-4-lh*0.9} textAnchor="start" fontSize={9.5*fontScale} fontWeight="700" fill="#111827" stroke="white" strokeWidth={2.6} style={{fontFamily:"var(--ff)",paintOrder:"stroke"}}>{r.name}</text>
              <text x={PLOT_L+3} y={bl-4} textAnchor="start" fontSize={8*fontScale} fontWeight="600" fill="#6B7280" stroke="white" strokeWidth={2.2} style={{fontFamily:"var(--ff)",paintOrder:"stroke"}}>{gmfiTxt}</text>
            </g>
          );
        })}
        {showPct&&pctPos==="peak"&&rows.map((r,i)=><text key={"pk"+i} x={r.labelX} y={r.labelY} textAnchor="start" fontSize={9*fontScale} fontWeight="700" fill={gateColor} stroke="white" strokeWidth={labelMode==="inline"?2.4:0} style={{fontFamily:"var(--ff)",paintOrder:"stroke"}}>{r.pct}%</text>)}
        <text x={PLOT_L+PLOT_W/2} y={lastRowBase+40} textAnchor="middle" fontSize="10.5" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{xLabel}</text>
      </svg>
    </div>
  );
}

// ─── Legend ──────────────────────────────────────────
const GripIcon=()=><svg width="12" height="16" viewBox="0 0 12 16" fill="currentColor" aria-hidden="true"><circle cx="3.5" cy="3" r="1.3"/><circle cx="8.5" cy="3" r="1.3"/><circle cx="3.5" cy="8" r="1.3"/><circle cx="8.5" cy="8" r="1.3"/><circle cx="3.5" cy="13" r="1.3"/><circle cx="8.5" cy="13" r="1.3"/></svg>;
const KebabIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor" aria-hidden="true"><circle cx="12" cy="5" r="1.7"/><circle cx="12" cy="12" r="1.7"/><circle cx="12" cy="19" r="1.7"/></svg>;
const EyeIcon=({off})=>off
  ?<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="M17.94 17.94A10.07 10.07 0 0112 20c-7 0-11-8-11-8a18.45 18.45 0 015.06-5.94M9.9 4.24A9.12 9.12 0 0112 4c7 0 11 8 11 8a18.5 18.5 0 01-2.16 3.19m-6.72-1.07a3 3 0 11-4.24-4.24"/><line x1="1" y1="1" x2="23" y2="23"/></svg>
  :<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"/><circle cx="12" cy="12" r="3"/></svg>;
const menuItemStyle={display:"block",width:"100%",textAlign:"left",padding:"7px 10px",border:"none",background:"none",borderRadius:6,fontSize:12.5,fontWeight:500,cursor:"pointer",fontFamily:"var(--ff)"};
function Legend({samples,colors,onColorChange,onReorder,hiddenFlags,onToggleHide,onDelete,onRename}){
  const[menu,setMenu]=useState(null);
  const[dragI,setDragI]=useState(null);
  const[overI,setOverI]=useState(null);
  const menuRef=useRef(null);
  useEffect(()=>{if(menu===null)return;const h=e=>{if(menuRef.current&&!menuRef.current.contains(e.target))setMenu(null);};document.addEventListener("mousedown",h);return()=>document.removeEventListener("mousedown",h);},[menu]);
  return <div style={{display:"flex",flexWrap:"wrap",gap:8,padding:"12px 16px",background:"white",borderRadius:10,border:"1px solid #E5E7EB"}}>
    <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em",lineHeight:"28px",marginRight:2}}>Legend</span>
    {samples.map((s,i)=>{const hid=hiddenFlags&&hiddenFlags[i];const isOver=overI===i&&dragI!==null&&dragI!==i;
      return <div key={i}
        onDragOver={onReorder?e=>{e.preventDefault();setOverI(i);}:undefined}
        onDragLeave={onReorder?()=>setOverI(o=>o===i?null:o):undefined}
        onDrop={onReorder?e=>{e.preventDefault();if(dragI!==null&&dragI!==i)onReorder(dragI,i);setDragI(null);setOverI(null);}:undefined}
        style={{display:"flex",alignItems:"center",gap:6,padding:"3px 6px 3px 3px",borderRadius:8,border:"1px solid "+(isOver?"#3B82F6":"#E5E7EB"),background:dragI===i?"#EFF6FF":"white",opacity:hid?0.45:1}}>
        {onReorder&&<span draggable title="Drag to reorder"
          onDragStart={e=>{setDragI(i);e.dataTransfer.effectAllowed="move";try{e.dataTransfer.setData("text/plain",String(i));}catch(_){}}}
          onDragEnd={()=>{setDragI(null);setOverI(null);}}
          style={{display:"flex",alignItems:"center",color:"#C4C9D0",cursor:"grab",padding:"0 1px"}}><GripIcon/></span>}
        <ColorPicker color={colors[i]||PALETTE[i%PALETTE.length]} onChange={c=>onColorChange(i,c)}/>
        {onRename
          ?<EditableText value={s.name} onChange={n=>onRename(i,n)} style={{fontSize:12,fontWeight:500,color:"#374151",textDecoration:hid?"line-through":"none"}} inputStyle={{fontSize:12,fontWeight:500,minWidth:80}}/>
          :<span style={{fontSize:12,color:"#374151",textDecoration:hid?"line-through":"none"}}>{s.name}</span>}
        {onToggleHide&&<button onClick={()=>onToggleHide(i)} title={hid?"Show on plots":"Hide from plots"}
          style={{display:"flex",alignItems:"center",justifyContent:"center",width:26,height:26,padding:0,borderRadius:6,border:"1px solid #E5E7EB",background:"white",cursor:"pointer",color:hid?"#3B82F6":"#6B7280"}}><EyeIcon off={hid}/></button>}
        {onDelete&&<div style={{position:"relative"}} ref={menu===i?menuRef:null}>
          <button onClick={()=>setMenu(menu===i?null:i)} title="More options"
            style={{display:"flex",alignItems:"center",justifyContent:"center",width:26,height:26,padding:0,borderRadius:6,border:"1px solid #E5E7EB",background:menu===i?"#F3F4F6":"white",cursor:"pointer",color:"#6B7280"}}><KebabIcon/></button>
          {menu===i&&<div style={{position:"absolute",top:"112%",right:0,zIndex:60,background:"white",border:"1px solid #E5E7EB",borderRadius:8,boxShadow:"0 8px 24px rgba(0,0,0,0.12)",padding:4,minWidth:130}}>
            {onToggleHide&&<button onClick={()=>{onToggleHide(i);setMenu(null);}} style={{...menuItemStyle,color:"#374151"}} onMouseEnter={e=>e.currentTarget.style.background="#F3F4F6"} onMouseLeave={e=>e.currentTarget.style.background="none"}>{hid?"Show on plots":"Hide from plots"}</button>}
            <button onClick={()=>{setMenu(null);if(window.confirm("Remove \""+s.name+"\" from the analysis?"))onDelete(i);}} style={{...menuItemStyle,color:"#EF4444"}} onMouseEnter={e=>e.currentTarget.style.background="#FEF2F2"} onMouseLeave={e=>e.currentTarget.style.background="none"}>Delete sample</button>
          </div>}
        </div>}
      </div>;})}
  </div>;
}

// ─── Shared dot-plot geometry + canvas painter ───────
// One source of truth so on-screen plots, single-plot export and whole-panel export all
// render identically (matching the reference FACS panel: rainbow dots, corner percentages).
const DPG={W:400,H:400,ML:60,MR:20,MT:36,MB:52};
function quadPct(xVals,yVals,gateX,gateY){
  let q1=0,q2=0,q3=0,q4=0,xp=0,yp=0;const n=Math.min(xVals.length,yVals.length);
  for(let i=0;i<n;i++){const xpos=xVals[i]>=gateX,ypos=yVals[i]>=gateY;if(xpos)xp++;if(ypos)yp++;if(!xpos&&ypos)q1++;else if(xpos&&ypos)q2++;else if(!xpos&&!ypos)q3++;else q4++;}
  const pc=v=>n?((v/n)*100).toFixed(1):"0.0";
  return{q1:pc(q1),q2:pc(q2),q3:pc(q3),q4:pc(q4),xPos:pc(xp),xNeg:pc(n-xp),yPos:pc(yp),yNeg:pc(n-yp),n};
}
// Paint a complete dot plot (dots, axes, ticks, titles+arrows, gate, corner %) into ctx at (ox,oy).
function paintDotPlot(ctx,ox,oy,o){
  const {W,H,ML,MR,MT,MB}=DPG,PWq=W-ML-MR,PHq=H-MT-MB;
  const {xVals,yVals,name,xCh,yCh,xRange,yRange,gateX,gateY,gateMode="quad",dotMode="color",color="#2563eb",
    dotSize=2,tickFont=1,labelFont=1,showGate=true,showPct=true,gateColor="#111827",xScale="log",yScale="log",
    showXAxis=true,showYAxis=true}=o;
  const isXLog=xScale==="log",isYLog=yScale==="log";
  const tfx=v=>isXLog?Math.log10(v>xRange.dMin?v:xRange.dMin):T(v);
  const tfy=v=>isYLog?Math.log10(v>yRange.dMin?v:yRange.dMin):T(v);
  const xS=v=>ox+ML+((tfx(v)-xRange.lo)/(xRange.hi-xRange.lo))*PWq;
  const yS=v=>oy+MT+PHq-((tfy(v)-yRange.lo)/(yRange.hi-yRange.lo))*PHq;
  const fmtX=isXLog?fmtLogTick:fmtTick,fmtY=isYLog?fmtLogTick:fmtTick;
  const xTicks=isXLog?logTicks(xRange.dMin,xRange.dMax,PWq):axisTicks(xRange.dMin,xRange.dMax,PWq);
  const yTicks=isYLog?logTicks(yRange.dMin,yRange.dMax,PHq):axisTicks(yRange.dMin,yRange.dMax,PHq);
  // dots (clipped to plot area)
  ctx.save();ctx.beginPath();ctx.rect(ox+ML,oy+MT,PWq,PHq);ctx.clip();
  const n=Math.min(xVals.length,yVals.length);
  const step=n>55000?Math.ceil(n/55000):1;
  const ds=Math.max(1,Math.round(dotSize)),off=(ds-1)>>1;
  if(dotMode==="density"){
    paintDensity(ctx,xVals,yVals,xS,yS,ox+ML,oy+MT,PWq,PHq,dotSize);
  }else{
    ctx.fillStyle=color;
    for(let i=0;i<n;i+=step)ctx.fillRect((xS(xVals[i])|0)-off,(yS(yVals[i])|0)-off,ds,ds);
  }
  ctx.restore();
  // axis lines (L)
  ctx.strokeStyle="#111111";ctx.lineWidth=1.1;
  ctx.beginPath();ctx.moveTo(ox+ML,oy+MT);ctx.lineTo(ox+ML,oy+MT+PHq);ctx.lineTo(ox+ML+PWq,oy+MT+PHq);ctx.stroke();
  // x ticks + labels
  ctx.fillStyle="#1F2937";ctx.font=(9.5*tickFont).toFixed(1)+"px sans-serif";ctx.textAlign="center";
  xTicks.forEach(({v,label})=>{const x=xS(v);ctx.strokeStyle=label?"#111111":"#6B7280";ctx.lineWidth=label?1.1:0.8;ctx.beginPath();ctx.moveTo(x,oy+MT+PHq);ctx.lineTo(x,oy+MT+PHq+(label?6:3.5));ctx.stroke();if(label&&showXAxis)ctx.fillText(fmtX(v),x,oy+MT+PHq+15+4*tickFont);});
  // y ticks + labels
  ctx.textAlign="right";
  yTicks.forEach(({v,label})=>{const y=yS(v);ctx.strokeStyle=label?"#111111":"#6B7280";ctx.lineWidth=label?1.1:0.8;ctx.beginPath();ctx.moveTo(ox+ML-(label?6:3.5),y);ctx.lineTo(ox+ML,y);ctx.stroke();if(label&&showYAxis)ctx.fillText(fmtY(v),ox+ML-9,y+4);});
  // axis titles
  if(showXAxis){
    ctx.fillStyle="#111827";ctx.font="600 "+(11*labelFont).toFixed(1)+"px sans-serif";ctx.textAlign="center";ctx.fillText(xCh,ox+ML+PWq/2,oy+H-14);
  }
  if(showYAxis){
    ctx.fillStyle="#111827";ctx.font="600 "+(11*labelFont).toFixed(1)+"px sans-serif";ctx.save();ctx.translate(ox+13,oy+MT+PHq/2);ctx.rotate(-Math.PI/2);ctx.textAlign="center";ctx.fillText(yCh,0,0);ctx.restore();
  }
  // gate lines
  const gx=xS(gateX),gy=yS(gateY);
  if(showGate){ctx.setLineDash([6,4]);ctx.strokeStyle=gateColor;ctx.lineWidth=1;if(gateMode!=="horiz"){ctx.beginPath();ctx.moveTo(gx,oy+MT);ctx.lineTo(gx,oy+MT+PHq);ctx.stroke();}if(gateMode!=="vert"){ctx.beginPath();ctx.moveTo(ox+ML,gy);ctx.lineTo(ox+ML+PWq,gy);ctx.stroke();}ctx.setLineDash([]);}
  // corner percentages (away from the dots)
  if(showPct){
    const q=quadPct(xVals,yVals,gateX,gateY);
    ctx.font="bold 13px sans-serif";ctx.fillStyle="#111827";
    const pad=9,Lx=ox+ML+pad,Rx=ox+ML+PWq-pad,Ty=oy+MT+pad+11,By=oy+MT+PHq-pad;
    if(gateMode==="quad"){
      ctx.textAlign="left";ctx.fillText(q.q1+"%",Lx,Ty);ctx.fillText(q.q3+"%",Lx,By);
      ctx.textAlign="right";ctx.fillText(q.q2+"%",Rx,Ty);ctx.fillText(q.q4+"%",Rx,By);
    }else if(gateMode==="vert"){
      ctx.textAlign="left";ctx.fillText(q.xNeg+"%",Lx,Ty);
      ctx.textAlign="right";ctx.fillText(q.xPos+"%",Rx,Ty);
    }else{
      ctx.textAlign="left";ctx.fillText(q.yPos+"%",Lx,Ty);ctx.fillText(q.yNeg+"%",Lx,By);
    }
  }
  // sample title
  ctx.font="bold 13px sans-serif";ctx.fillStyle="#111827";ctx.textAlign="center";ctx.fillText(name,ox+W/2,oy+18);
}

// ─── Dot Plot Component (Canvas) ─────────────────────
function DotPlot({xVals,yVals,name,xCh,yCh,xRange,yRange,gateX,gateY,onGateChange,color,onNameChange,gateColor,showGate=true,showPct=true,gateMode="quad",dotMode="color",xScale="log",yScale="log",tickFont=1,labelFont=1,dotSize=2,showStats=true,bare=false}){
  const canvasRef=useRef(null);
  const containerRef=useRef(null);
  const dragRef=useRef(null);

  const W=400,H=400;
  const ML=60,MR=20,MT=36,MB=52;
  const PWq=W-ML-MR;
  const PHq=H-MT-MB;

  const isXLog=xScale==="log",isYLog=yScale==="log";
  const tfx=useCallback(v=>isXLog?Math.log10(v>xRange.dMin?v:xRange.dMin):T(v),[isXLog,xRange]);
  const tfy=useCallback(v=>isYLog?Math.log10(v>yRange.dMin?v:yRange.dMin):T(v),[isYLog,yRange]);
  const xS=useCallback(v=>ML+((tfx(v)-xRange.lo)/(xRange.hi-xRange.lo))*PWq,[tfx,xRange]);
  const yS=useCallback(v=>MT+PHq-((tfy(v)-yRange.lo)/(yRange.hi-yRange.lo))*PHq,[tfy,yRange]);
  const fmtX=isXLog?fmtLogTick:fmtTick,fmtY=isYLog?fmtLogTick:fmtTick;

  const xTicks=useMemo(()=>isXLog?logTicks(xRange.dMin,xRange.dMax,PWq):axisTicks(xRange.dMin,xRange.dMax,PWq),[isXLog,xRange]);
  const yTicks=useMemo(()=>isYLog?logTicks(yRange.dMin,yRange.dMax,PHq):axisTicks(yRange.dMin,yRange.dMax,PHq),[isYLog,yRange]);

  const quadrants=useMemo(()=>{
    let q1=0,q2=0,q3=0,q4=0,xp=0,yp=0;
    const n=Math.min(xVals.length,yVals.length);
    for(let i=0;i<n;i++){
      const xpos=xVals[i]>=gateX,ypos=yVals[i]>=gateY;
      if(xpos)xp++;if(ypos)yp++;
      if(!xpos&&ypos)q1++;
      else if(xpos&&ypos)q2++;
      else if(!xpos&&!ypos)q3++;
      else q4++;
    }
    const pc=v=>n?((v/n)*100).toFixed(1):"0.0";
    return{q1:pc(q1),q2:pc(q2),q3:pc(q3),q4:pc(q4),xPos:pc(xp),xNeg:pc(n-xp),yPos:pc(yp),yNeg:pc(n-yp),n};
  },[xVals,yVals,gateX,gateY]);

  useEffect(()=>{
    const canvas=canvasRef.current;
    if(!canvas)return;
    const ctx=canvas.getContext("2d");
    const dpr=window.devicePixelRatio||1;
    canvas.width=W*dpr;
    canvas.height=H*dpr;
    ctx.scale(dpr,dpr);
    ctx.clearRect(0,0,W,H);
    ctx.save();
    ctx.beginPath();
    ctx.rect(ML,MT,PWq,PHq);
    ctx.clip();
    const n=Math.min(xVals.length,yVals.length);
    const maxPts=55000;
    const step=n>maxPts?Math.ceil(n/maxPts):1;
    const ds=Math.max(1,Math.round(dotSize)),off=(ds-1)>>1;
    if(dotMode==="density"){
      ctx.globalAlpha=1;
      paintDensity(ctx,xVals,yVals,xS,yS,ML,MT,PWq,PHq,dotSize);
    }else{
      ctx.globalAlpha=1;
      ctx.fillStyle=color;
      for(let i=0;i<n;i+=step){ctx.fillRect((xS(xVals[i])|0)-off,(yS(yVals[i])|0)-off,ds,ds);}
    }
    ctx.restore();
  },[xVals,yVals,xRange,yRange,color,xS,yS,dotMode,dotSize]);

  const gxPx=xS(gateX);
  const gyPx=yS(gateY);

  const getValFromMouse=useCallback((e,axis)=>{
    const rect=containerRef.current.getBoundingClientRect();
    if(axis==="x"){
      const px=(e.clientX-rect.left)/rect.width*W;
      const frac=(px-ML)/PWq;
      const lv=xRange.lo+Math.max(0,Math.min(1,frac))*(xRange.hi-xRange.lo);
      return Math.round(isXLog?10**lv:invT(lv));
    }else{
      const py=(e.clientY-rect.top)/rect.height*H;
      const frac=1-(py-MT)/PHq;
      const lv=yRange.lo+Math.max(0,Math.min(1,frac))*(yRange.hi-yRange.lo);
      return Math.round(isYLog?10**lv:invT(lv));
    }
  },[xRange,yRange,isXLog,isYLog]);

  const onMouseDown=useCallback(e=>{
    const rect=containerRef.current.getBoundingClientRect();
    const mx=(e.clientX-rect.left)/rect.width*W;
    const my=(e.clientY-rect.top)/rect.height*H;
    const nearX=gateMode!=="horiz"&&Math.abs(mx-gxPx)<10;
    const nearY=gateMode!=="vert"&&Math.abs(my-gyPx)<10;
    if(gateMode==="quad"&&nearX&&nearY)dragRef.current="both";
    else if(nearX)dragRef.current="x";
    else if(nearY)dragRef.current="y";
    else dragRef.current=null;
    if(dragRef.current)e.preventDefault();
  },[gxPx,gyPx,gateMode]);

  const onMouseMove=useCallback(e=>{
    if(!dragRef.current)return;
    const d=dragRef.current;
    const newGateX=(d==="x"||d==="both")?getValFromMouse(e,"x"):gateX;
    const newGateY=(d==="y"||d==="both")?getValFromMouse(e,"y"):gateY;
    onGateChange(newGateX,newGateY);
  },[getValFromMouse,gateX,gateY,onGateChange]);

  const onMouseUp=useCallback(()=>{dragRef.current=null;},[]);

  const[cursor,setCursor]=useState("default");
  const onHover=useCallback(e=>{
    if(dragRef.current)return;
    const rect=containerRef.current.getBoundingClientRect();
    const mx=(e.clientX-rect.left)/rect.width*W;
    const my=(e.clientY-rect.top)/rect.height*H;
    const nearX=gateMode!=="horiz"&&Math.abs(mx-gxPx)<10;
    const nearY=gateMode!=="vert"&&Math.abs(my-gyPx)<10;
    if(gateMode==="quad"&&nearX&&nearY)setCursor("move");
    else if(nearX)setCursor("ew-resize");
    else if(nearY)setCursor("ns-resize");
    else setCursor("default");
  },[gxPx,gyPx,gateMode]);

  const exportPlot=useCallback(()=>{
    const scale=3;
    const expCanvas=document.createElement("canvas");
    expCanvas.width=W*scale;expCanvas.height=H*scale;
    const ctx=expCanvas.getContext("2d");
    ctx.scale(scale,scale);
    ctx.fillStyle="white";ctx.fillRect(0,0,W,H);
    paintDotPlot(ctx,0,0,{xVals,yVals,name,xCh,yCh,xRange,yRange,gateX,gateY,gateMode,dotMode,color,dotSize,tickFont,labelFont,showGate,showPct,gateColor,xScale,yScale,showXAxis:true,showYAxis:true});
    downloadDataURL(expCanvas.toDataURL("image/png"),name.replace(/\s+/g,"_")+".png");
  },[xVals,yVals,gateX,gateY,name,xCh,yCh,xRange,yRange,color,dotMode,dotSize,showGate,showPct,gateColor,gateMode,xScale,yScale,tickFont,labelFont]);

  const cornerStyle={fontFamily:"var(--ff)",paintOrder:"stroke"};
  const pad=9,Lx=ML+pad,Rx=ML+PWq-pad,Ty=MT+pad+11,By=MT+PHq-pad;
  return(
    <div style={bare?{background:"transparent"}:{background:"white",borderRadius:10,border:"1px solid #e2e5ea",padding:"10px 6px 6px",boxShadow:"0 1px 4px rgba(0,0,0,0.04)"}}>
      <div style={{display:"flex",alignItems:"center",justifyContent:"center",gap:8,marginBottom:4}}>
        {onNameChange
          ?<EditableText value={name} onChange={onNameChange} style={{fontSize:13,fontWeight:700,color:"#111827"}} inputStyle={{fontSize:13,fontWeight:700,textAlign:"center",width:"70%"}}/>
          :<span style={{fontSize:13,fontWeight:700,color:"#111827",fontFamily:"var(--ff)"}}>{name}</span>}
        <button onClick={exportPlot} title="Export this plot as PNG"
          style={{background:"none",border:"none",cursor:"pointer",padding:2,color:"#CBD5E1",display:"flex"}}
          onMouseEnter={e=>{e.currentTarget.style.color="#6B7280";}}
          onMouseLeave={e=>{e.currentTarget.style.color="#CBD5E1";}}>
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
        </button>
      </div>
      <div ref={containerRef} style={{position:"relative",width:W,maxWidth:"100%",aspectRatio:W+"/"+H,margin:"0 auto"}}
        onMouseDown={onMouseDown} onMouseMove={e=>{onMouseMove(e);onHover(e);}} onMouseUp={onMouseUp} onMouseLeave={onMouseUp}>
        <canvas ref={canvasRef} style={{position:"absolute",top:0,left:0,width:"100%",height:"100%",cursor}}/>
        <svg viewBox={"0 0 "+W+" "+H} style={{position:"absolute",top:0,left:0,width:"100%",height:"100%",pointerEvents:"none"}}>
          <line x1={ML} x2={ML} y1={MT} y2={MT+PHq} stroke="#111111" strokeWidth={1.1}/>
          <line x1={ML} x2={ML+PWq} y1={MT+PHq} y2={MT+PHq} stroke="#111111" strokeWidth={1.1}/>
          {xTicks.map(({v,label})=>{const x=xS(v);return(<g key={"x"+v}><line x1={x} x2={x} y1={MT+PHq} y2={MT+PHq+(label?6:3.5)} stroke={label?"#111111":"#6B7280"} strokeWidth={label?1.1:0.8}/>{label&&<text x={x} y={MT+PHq+15+4*tickFont} textAnchor="middle" fontSize={9.5*tickFont} fill="#1F2937" style={{fontFamily:"var(--ff)"}}>{fmtX(v)}</text>}</g>);})}
          <text x={ML+PWq/2} y={H-14} textAnchor="middle" fontSize={11*labelFont} fill="#111827" fontWeight="600" style={{fontFamily:"var(--ff)"}}>{xCh}</text>
          {yTicks.map(({v,label})=>{const y=yS(v);return(<g key={"y"+v}><line x1={ML-(label?6:3.5)} x2={ML} y1={y} y2={y} stroke={label?"#111111":"#6B7280"} strokeWidth={label?1.1:0.8}/>{label&&<text x={ML-9} y={y+3.5} textAnchor="end" fontSize={9.5*tickFont} fill="#1F2937" style={{fontFamily:"var(--ff)"}}>{fmtY(v)}</text>}</g>);})}
          <text transform={"translate(13,"+(MT+PHq/2)+") rotate(-90)"} textAnchor="middle" fontSize={11*labelFont} fill="#111827" fontWeight="600" style={{fontFamily:"var(--ff)"}}>{yCh}</text>
          {showGate&&gateMode!=="horiz"&&<line x1={gxPx} x2={gxPx} y1={MT} y2={MT+PHq} stroke={gateColor} strokeWidth={1} strokeDasharray="6,4"/>}
          {showGate&&gateMode!=="vert"&&<line x1={ML} x2={ML+PWq} y1={gyPx} y2={gyPx} stroke={gateColor} strokeWidth={1} strokeDasharray="6,4"/>}
          {showPct&&gateMode==="quad"&&<g fontWeight="700" fontSize="13" fill="#111827" stroke="white" strokeWidth={3} style={cornerStyle}>
            <text x={Lx} y={Ty} textAnchor="start">{quadrants.q1}%</text>
            <text x={Rx} y={Ty} textAnchor="end">{quadrants.q2}%</text>
            <text x={Lx} y={By} textAnchor="start">{quadrants.q3}%</text>
            <text x={Rx} y={By} textAnchor="end">{quadrants.q4}%</text>
          </g>}
          {showPct&&gateMode==="vert"&&<g fontWeight="700" fontSize="13" fill="#111827" stroke="white" strokeWidth={3} style={cornerStyle}>
            <text x={Lx} y={Ty} textAnchor="start">{quadrants.xNeg}%</text>
            <text x={Rx} y={Ty} textAnchor="end">{quadrants.xPos}%</text>
          </g>}
          {showPct&&gateMode==="horiz"&&<g fontWeight="700" fontSize="13" fill="#111827" stroke="white" strokeWidth={3} style={cornerStyle}>
            <text x={Lx} y={Ty} textAnchor="start">{quadrants.yPos}%</text>
            <text x={Lx} y={By} textAnchor="start">{quadrants.yNeg}%</text>
          </g>}
        </svg>
      </div>
    </div>
  );
}

// ─── Icon components ─────────────────────────────────
const GridIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><rect x="3" y="3" width="7" height="7"/><rect x="14" y="3" width="7" height="7"/><rect x="3" y="14" width="7" height="7"/><rect x="14" y="14" width="7" height="7"/></svg>;
const RidgeIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><path d="M2 6h20"/><path d="M2 12h20"/><path d="M2 18h20"/><path d="M6 4l4 4-3 2 5 2"/><path d="M6 10l4 4-3 2 5 2"/></svg>;
const OverlayIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="M3 18c3-1 4-9 7-9s4 8 7 9"/><path d="M6 18c3-1 4-6 7-6s5 5 8 6"/></svg>;
const HistIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round"><path d="M3 20h18"/><path d="M5 20V10l3 4 3-8 3 6 3-10 3 8"/></svg>;
const QuadIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round"><circle cx="6" cy="6" r="2"/><circle cx="18" cy="18" r="2"/><circle cx="14" cy="8" r="2"/><circle cx="8" cy="16" r="2"/><circle cx="17" cy="11" r="2"/></svg>;

// ─── Shared Upload Box ───────────────────────────────
function UploadBox({samples,dragOver,setDragOver,onDrop,fileRef,handleFiles,accent,hint}){
  return(
    <div onDragOver={e=>{e.preventDefault();setDragOver(true);}} onDragLeave={()=>setDragOver(false)} onDrop={onDrop}
      onClick={()=>fileRef.current&&fileRef.current.click()}
      style={{flex:"1 1 260px",minHeight:96,borderRadius:12,border:"2px dashed "+(dragOver?accent:"#D1D5DB"),background:dragOver?accent+"14":"white",display:"flex",flexDirection:"column",alignItems:"center",justifyContent:"center",cursor:"pointer",transition:"all 0.15s",padding:"14px 20px"}}>
      <svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke={dragOver?accent:"#9CA3AF"} strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg>
      <span style={{fontSize:13,fontWeight:600,color:dragOver?accent:"#374151",marginTop:5}}>{samples.length>0?"+ Add more CSV files":"Drop CSV files or click to browse"}</span>
      <span style={{fontSize:11,color:"#9CA3AF",marginTop:3,textAlign:"center",lineHeight:1.5}}>{hint||<>One CSV per sample, <b>pre-gated to your population</b> (e.g. live singlets).<br/>FlowJo → Export → CSV (Channel Values), including each channel you'll plot.</>}</span>
      <input ref={fileRef} type="file" accept=".csv,.tsv,.txt" multiple onChange={e=>{handleFiles([...e.target.files]);e.target.value="";}} style={{display:"none"}}/>
    </div>
  );
}

const labelStyle={fontSize:10,fontWeight:600,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.05em"};
const inputStyle={marginTop:3,padding:"6px 8px",borderRadius:6,border:"1px solid #D1D5DB",fontSize:12.5,fontFamily:"var(--ff)",color:"#111827",boxSizing:"border-box"};

// ─── Histogram Mode ──────────────────────────────────
function HistogramMode({samples,allHeaders,colors,updateSampleName,updateColor,removeSample,swapSamples,moveSample,hidden,toggleHidden,uploadProps,gateColor,showGate,onToggleGate,showPct,onTogglePct,initState,onState}){
  const S=initState||{};
  const[channel,setChannel]=useState(S.channel||"");
  const[gate,setGate]=useState(S.gate??500);
  const[gateLabel,setGateLabel]=useState(S.gateLabel??"PE+");
  const[xLabel,setXLabel]=useState(S.xLabel||"");
  const[yLabel,setYLabel]=useState(S.yLabel??"Count");
  const[xMinInput,setXMinInput]=useState(S.xMinInput??"");
  const[xMaxInput,setXMaxInput]=useState(S.xMaxInput??"");
  const[appliedXMin,setAppliedXMin]=useState(S.appliedXMin??"");
  const[appliedXMax,setAppliedXMax]=useState(S.appliedXMax??"");
  const[viewMode,setViewMode]=useState(S.viewMode??"grid"); // grid | ridge | overlay
  const[overlaySel,setOverlaySel]=useState(S.overlaySel!==undefined?S.overlaySel:null); // null = all
  const[overlayNorm,setOverlayNorm]=useState(S.overlayNorm??true);
  const[yMode,setYMode]=useState(S.yMode??"count"); // "count" | "pct"
  const[log2Basis,setLog2Basis]=useState(S.log2Basis??"all"); // "all" | "gated"
  const[log2RefIdx,setLog2RefIdx]=useState(S.log2RefIdx!==undefined?S.log2RefIdx:null); // reference sample._idx; null = first
  const[ridgePctPos,setRidgePctPos]=useState(S.ridgePctPos??"left");
  const[overlayPctPos,setOverlayPctPos]=useState(S.overlayPctPos??"left");

  // Report state up for session saving
  useEffect(()=>{
    if(onState)onState({channel,gate,gateLabel,xLabel,yLabel,xMinInput,xMaxInput,appliedXMin,appliedXMax,viewMode,overlaySel,overlayNorm,yMode,log2Basis,log2RefIdx,ridgePctPos,overlayPctPos});
  });

  // Auto-pick channel once headers exist (or re-pick if current one was merged away)
  useEffect(()=>{
    if((!channel||!allHeaders.includes(channel))&&allHeaders.length){const det=detectPE(allHeaders);setChannel(det);setXLabel(det);}
  },[allHeaders,channel]);

  const hiddenSet=hidden||[];
  const allValid=useMemo(()=>samples.map((s,i)=>({...s,_idx:i})).filter(s=>s.columns[channel]&&s.columns[channel].length>0),[samples,channel]);
  const validSamples=useMemo(()=>allValid.filter(s=>!hiddenSet.includes(s._idx)),[allValid,hidden]);
  const gridCols=validSamples.length===1?1:validSamples.length<=4?2:3;
  const handleChannelChange=ch=>{setChannel(ch);setXLabel(ch);};
  const autoXAxisRange=useMemo(()=>analyzePooledValues(validSamples,channel),[validSamples,channel]);
  const parsedXAxisDraft=useMemo(()=>{
    const minRaw=xMinInput.trim();const maxRaw=xMaxInput.trim();
    if(minRaw===""&&maxRaw==="")return{kind:"auto"};
    if(minRaw===""||maxRaw==="")return{kind:"invalid"};
    const min=Number(minRaw);const max=Number(maxRaw);
    if(!Number.isFinite(min)||!Number.isFinite(max)||max<=min)return{kind:"invalid"};
    return{kind:"manual",min,max};
  },[xMinInput,xMaxInput]);
  const xDomain=useMemo(()=>{
    const min=Number(appliedXMin);const max=Number(appliedXMax);
    if(appliedXMin.trim()===""||appliedXMax.trim()==="")return null;
    if(!Number.isFinite(min)||!Number.isFinite(max)||max<=min)return null;
    return{lo:T(min),hi:T(max),dMin:min,dMax:max};
  },[appliedXMin,appliedXMax]);
  const axisModeLabel=xDomain?"Manual":"Auto";
  const hasPendingXAxis=xMinInput!==appliedXMin||xMaxInput!==appliedXMax;
  const applyXAxis=()=>{
    if(parsedXAxisDraft.kind==="auto"){setAppliedXMin("");setAppliedXMax("");return;}
    if(parsedXAxisDraft.kind!=="manual")return;
    setAppliedXMin(xMinInput.trim());setAppliedXMax(xMaxInput.trim());
  };
  const resetXAxis=()=>{setXMinInput("");setXMaxInput("");setAppliedXMin("");setAppliedXMax("");};
  const handleXAxisKeyDown=e=>{if(e.key==="Enter"){e.preventDefault();applyXAxis();}};

  // Overlay selection ─ effective set of sample _idx that are shown
  const selectedIdx=useMemo(()=>{
    if(overlaySel===null)return new Set(validSamples.map(s=>s._idx));
    return new Set(validSamples.filter(s=>overlaySel.includes(s._idx)).map(s=>s._idx));
  },[overlaySel,validSamples]);
  const overlaySamples=useMemo(()=>validSamples.filter(s=>selectedIdx.has(s._idx)),[validSamples,selectedIdx]);
  const toggleOverlay=idx=>{
    const base=overlaySel===null?validSamples.map(s=>s._idx):overlaySel;
    const next=base.includes(idx)?base.filter(i=>i!==idx):[...base,idx];
    setOverlaySel(next);
  };
  const overlayAll=()=>setOverlaySel(null);
  const overlayNone=()=>setOverlaySel([]);

  return(
    <>
      {/* Upload + Controls */}
      <div style={{maxWidth:1200,margin:"0 auto 20px",display:"flex",gap:14,flexWrap:"wrap",alignItems:"stretch"}}>
        <UploadBox {...uploadProps} samples={samples} accent="#3B82F6"/>
        {samples.length>0&&(
          <div style={{flex:"0 0 340px",background:"white",borderRadius:12,border:"1px solid #E5E7EB",padding:"14px 18px",display:"flex",flexDirection:"column",gap:12}}>
            <div style={{display:"flex",gap:12}}>
              <div style={{flex:1}}>
                <label style={labelStyle}>Channel</label>
                <select value={channel} onChange={e=>handleChannelChange(e.target.value)} style={{...inputStyle,width:"100%",background:"white"}}>
                  {allHeaders.filter(h=>!/Time/i.test(h)).map(h=><option key={h} value={h}>{h}</option>)}
                </select>
              </div>
              <div>
                <label style={labelStyle}>Gate</label>
                <input type="number" value={gate} onChange={e=>setGate(Number(e.target.value))} style={{...inputStyle,width:80}}/>
              </div>
              <div>
                <label style={labelStyle}>Label</label>
                <input value={gateLabel} onChange={e=>setGateLabel(e.target.value)} style={{...inputStyle,width:70}} placeholder="PE+"/>
              </div>
            </div>
            <div style={{display:"flex",gap:12}}>
              <div style={{flex:1}}>
                <label style={labelStyle}>X-axis</label>
                <input value={xLabel} onChange={e=>setXLabel(e.target.value)} style={{...inputStyle,width:"100%"}}/>
              </div>
              <div style={{flex:1}}>
                <label style={labelStyle}>Y-axis</label>
                <input value={yLabel} onChange={e=>setYLabel(e.target.value)} disabled={yMode==="pct"} placeholder={yMode==="pct"?"% of max":"Count"}
                  style={{...inputStyle,width:"100%",background:yMode==="pct"?"#F9FAFB":"white",color:yMode==="pct"?"#9CA3AF":"#111827"}}/>
              </div>
            </div>
            <div style={{display:"flex",alignItems:"center",gap:8}}>
              <label style={labelStyle}>Y scale</label>
              <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
                <button onClick={()=>setYMode("count")} style={{padding:"4px 12px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:yMode==="count"?"#EFF6FF":"white",color:yMode==="count"?"#3B82F6":"#9CA3AF"}} title="Raw event counts">Count</button>
                <button onClick={()=>setYMode("pct")} style={{padding:"4px 12px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:yMode==="pct"?"#EFF6FF":"white",color:yMode==="pct"?"#3B82F6":"#9CA3AF"}} title="Each histogram scaled to its own peak (0–100%)">% max</button>
              </div>
              <span style={{fontSize:10.5,color:"#9CA3AF"}}>Grid histograms</span>
            </div>
            <div>
              <div style={{display:"flex",alignItems:"center",justifyContent:"space-between"}}>
                <label style={labelStyle}>X-range</label>
                <div style={{display:"flex",alignItems:"center",gap:8}}>
                  <span style={{fontSize:10.5,fontWeight:600,color:xDomain?"#2563EB":"#9CA3AF"}}>{axisModeLabel}</span>
                  <button onClick={resetXAxis} style={{padding:"3px 8px",borderRadius:999,border:"1px solid #E5E7EB",background:"white",fontSize:10.5,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>Auto</button>
                </div>
              </div>
              <div style={{display:"flex",gap:12,marginTop:3}}>
                <input value={xMinInput} onChange={e=>setXMinInput(e.target.value)} onKeyDown={handleXAxisKeyDown} placeholder={String(Math.round(autoXAxisRange.dMin))} style={{...inputStyle,width:"100%"}}/>
                <input value={xMaxInput} onChange={e=>setXMaxInput(e.target.value)} onKeyDown={handleXAxisKeyDown} placeholder={String(Math.round(autoXAxisRange.dMax))} style={{...inputStyle,width:"100%"}}/>
              </div>
              <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",marginTop:4,fontSize:10.5,color:"#9CA3AF",gap:12}}>
                <span>Min / max in raw channel units</span>
                <div style={{display:"flex",alignItems:"center",gap:8}}>
                  {parsedXAxisDraft.kind==="invalid"&&<span style={{color:"#DC2626"}}>Enter both values with max greater than min</span>}
                  {!xDomain&&parsedXAxisDraft.kind!=="invalid"&&<span>{"auto "+Math.round(autoXAxisRange.dMin)+" to "+Math.round(autoXAxisRange.dMax)}</span>}
                  <button onClick={applyXAxis} disabled={parsedXAxisDraft.kind==="invalid"||!hasPendingXAxis}
                    style={{padding:"4px 10px",borderRadius:6,border:"1px solid "+(parsedXAxisDraft.kind==="invalid"||!hasPendingXAxis?"#E5E7EB":"#BFDBFE"),background:parsedXAxisDraft.kind==="invalid"||!hasPendingXAxis?"#F9FAFB":"#EFF6FF",fontSize:11,fontWeight:700,color:parsedXAxisDraft.kind==="invalid"||!hasPendingXAxis?"#9CA3AF":"#2563EB",cursor:parsedXAxisDraft.kind==="invalid"||!hasPendingXAxis?"default":"pointer",fontFamily:"var(--ff)"}}>Apply</button>
                </div>
              </div>
            </div>
            <div style={{display:"flex",alignItems:"center",gap:8,marginTop:"auto"}}>
              <span style={{fontSize:12,color:"#6B7280"}}>{validSamples.length+" sample"+(validSamples.length!==1?"s":"")}</span>
              {validSamples.length>1&&(
                <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden",marginLeft:4}}>
                  <button onClick={()=>setViewMode("grid")} style={{padding:"3px 9px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",display:"flex",alignItems:"center",gap:4,background:viewMode==="grid"?"#EFF6FF":"white",color:viewMode==="grid"?"#3B82F6":"#9CA3AF"}}><GridIcon/> Grid</button>
                  <button onClick={()=>setViewMode("overlay")} style={{padding:"3px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",display:"flex",alignItems:"center",gap:4,background:viewMode==="overlay"?"#EFF6FF":"white",color:viewMode==="overlay"?"#3B82F6":"#9CA3AF"}}><OverlayIcon/> Overlay</button>
                  <button onClick={()=>setViewMode("ridge")} style={{padding:"3px 9px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",display:"flex",alignItems:"center",gap:4,background:viewMode==="ridge"?"#EFF6FF":"white",color:viewMode==="ridge"?"#3B82F6":"#9CA3AF"}}><RidgeIcon/> Ridge</button>
                </div>
              )}
            </div>
          </div>
        )}
      </div>

      {/* Legend */}
      {allValid.length>1&&(
        <div style={{maxWidth:1200,margin:"0 auto 16px"}}>
          <Legend samples={allValid} colors={allValid.map(s=>colors[s._idx]||PALETTE[s._idx%PALETTE.length])} onColorChange={(vi,c)=>updateColor(allValid[vi]._idx,c)}
            onReorder={moveSample?(fromVi,toVi)=>moveSample(allValid[fromVi]._idx,allValid[toVi]._idx):undefined}
            onRename={(vi,n)=>updateSampleName(allValid[vi]._idx,n)}
            hiddenFlags={allValid.map(s=>hiddenSet.includes(s._idx))}
            onToggleHide={vi=>toggleHidden(allValid[vi]._idx)}
            onDelete={vi=>removeSample(allValid[vi]._idx)}/>
        </div>
      )}

      {/* Grid View */}
      {validSamples.length>0&&viewMode==="grid"&&(
        <div style={{maxWidth:1200,margin:"0 auto",display:"grid",gridTemplateColumns:"repeat("+gridCols+", 1fr)",gap:16}}>
          {validSamples.map(s=>(
            <div key={s._idx} style={{position:"relative"}}>
              <button onClick={()=>removeSample(s._idx)} title="Remove sample"
                style={{position:"absolute",top:4,right:4,zIndex:10,width:22,height:22,borderRadius:6,border:"none",background:"transparent",color:"#CBD5E1",cursor:"pointer",fontSize:16,lineHeight:"20px",fontFamily:"var(--ff)",display:"flex",alignItems:"center",justifyContent:"center"}}
                onMouseEnter={e=>{e.currentTarget.style.background="#FEE2E2";e.currentTarget.style.color="#EF4444";}} onMouseLeave={e=>{e.currentTarget.style.background="transparent";e.currentTarget.style.color="#CBD5E1";}}>×</button>
              <Histogram values={s.columns[channel]} name={s.name} color={colors[s._idx]||PALETTE[s._idx%PALETTE.length]}
                xLabel={xLabel} yLabel={yLabel} gateValue={gate} onGateChange={setGate} onNameChange={n=>updateSampleName(s._idx,n)} xDomain={xDomain||autoXAxisRange} gateLabel={gateLabel} gateColor={gateColor} yMode={yMode} showGate={showGate} showPct={showPct}/>
            </div>
          ))}
        </div>
      )}

      {/* Overlay View */}
      {validSamples.length>1&&viewMode==="overlay"&&(
        <div style={{maxWidth:1200,margin:"0 auto"}}>
          {/* Selector bar */}
          <div style={{background:"white",borderRadius:10,border:"1px solid #E5E7EB",padding:"10px 14px",marginBottom:14,display:"flex",flexWrap:"wrap",alignItems:"center",gap:"8px 10px"}}>
            <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em",marginRight:2}}>Overlay</span>
            {validSamples.map(s=>{
              const on=selectedIdx.has(s._idx);
              const c=colors[s._idx]||PALETTE[s._idx%PALETTE.length];
              return(
                <button key={s._idx} onClick={()=>toggleOverlay(s._idx)}
                  style={{display:"flex",alignItems:"center",gap:6,padding:"4px 10px",borderRadius:999,cursor:"pointer",fontFamily:"var(--ff)",fontSize:12,fontWeight:600,
                    border:"1px solid "+(on?c:"#E5E7EB"),background:on?c+"14":"white",color:on?"#111827":"#9CA3AF",transition:"all 0.12s"}}>
                  <span style={{width:11,height:11,borderRadius:3,background:on?c:"#D1D5DB",flexShrink:0}}/>
                  {s.name}
                </button>
              );
            })}
            <div style={{marginLeft:"auto",display:"flex",gap:6,alignItems:"center"}}>
              <button onClick={overlayAll} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>All</button>
              <button onClick={overlayNone} style={{padding:"4px 10px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:11,fontWeight:600,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>None</button>
              <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
                <button onClick={()=>setOverlayNorm(true)} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:overlayNorm?"#EFF6FF":"white",color:overlayNorm?"#3B82F6":"#9CA3AF"}} title="Scale each curve to its own peak">Modal</button>
                <button onClick={()=>setOverlayNorm(false)} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:!overlayNorm?"#EFF6FF":"white",color:!overlayNorm?"#3B82F6":"#9CA3AF"}} title="Shared count axis across curves">Counts</button>
              </div>
            </div>
          </div>
          {overlaySamples.length>0
            ?<OverlayHistogram samples={overlaySamples} colors={overlaySamples.map(s=>colors[s._idx]||PALETTE[s._idx%PALETTE.length])}
              channel={channel} xLabel={xLabel} yLabel={yLabel} gateValue={gate} onGateChange={setGate} xDomain={xDomain} gateLabel={gateLabel} normalize={overlayNorm} gateColor={gateColor} showGate={showGate} showPct={showPct} onToggleGate={onToggleGate} onTogglePct={onTogglePct} pctPos={overlayPctPos} setPctPos={setOverlayPctPos}/>
            :<div style={{textAlign:"center",color:"#9CA3AF",fontSize:13,padding:"40px 0",background:"white",borderRadius:12,border:"1px solid #E5E7EB"}}>Select at least one sample above to overlay.</div>}
        </div>
      )}

      {/* Ridge View */}
      {validSamples.length>1&&viewMode==="ridge"&&(
        <div style={{maxWidth:1000,margin:"0 auto"}}>
          <RidgePlot samples={validSamples} colors={validSamples.map(s=>colors[s._idx]||PALETTE[s._idx%PALETTE.length])}
            channel={channel} gateValue={gate} onGateChange={setGate} xLabel={xLabel} xDomain={xDomain} gateLabel={gateLabel} gateColor={gateColor} showGate={showGate} onToggleGate={onToggleGate} showPct={showPct} onTogglePct={onTogglePct} pctPos={ridgePctPos} setPctPos={setRidgePctPos}/>
        </div>
      )}

      {/* Summary table + log2 fold-change */}
      {validSamples.length>1&&(()=>{
        const rows=validSamples.map(s=>{
          const vals=s.columns[channel];
          let pc=0,ps=0,gLogSum=0,gLogN=0,allLogSum=0,allLogN=0;
          for(let k=0;k<vals.length;k++){
            if(vals[k]>0){allLogSum+=Math.log(vals[k]);allLogN++;}
            if(vals[k]>=gate){pc++;ps+=vals[k];if(vals[k]>0){gLogSum+=Math.log(vals[k]);gLogN++;}}
          }
          return{s,n:vals.length,pct:((pc/vals.length)*100).toFixed(1),
            mfiNum:pc>0?ps/pc:null,
            gAll:allLogN>0?Math.exp(allLogSum/allLogN):null,
            gGated:gLogN>0?Math.exp(gLogSum/gLogN):null};
        });
        const ref=rows.find(r=>r.s._idx===log2RefIdx)||rows[0];
        const refG=log2Basis==="all"?ref.gAll:ref.gGated;
        const basisLabel=log2Basis==="all"?"all":gateLabel;
        return(
        <div style={{maxWidth:1200,margin:"20px auto 0"}}>
          <div style={{display:"flex",alignItems:"center",gap:10,flexWrap:"wrap",marginBottom:8,padding:"0 2px"}}>
            <span style={{fontSize:11.5,fontWeight:700,color:"#374151"}}>log₂ fold-change of gMFI vs</span>
            <select value={ref.s._idx} onChange={e=>setLog2RefIdx(Number(e.target.value))} style={{...inputStyle,marginTop:0,padding:"4px 8px",background:"white",fontSize:12}}>
              {rows.map(r=><option key={r.s._idx} value={r.s._idx}>{r.s.name}</option>)}
            </select>
            <span style={{fontSize:11,color:"#9CA3AF"}}>using</span>
            <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
              <button onClick={()=>setLog2Basis("all")} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:log2Basis==="all"?"#EFF6FF":"white",color:log2Basis==="all"?"#3B82F6":"#9CA3AF"}}>gMFI (all)</button>
              <button onClick={()=>setLog2Basis("gated")} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:log2Basis==="gated"?"#EFF6FF":"white",color:log2Basis==="gated"?"#3B82F6":"#9CA3AF"}}>gMFI ({gateLabel})</button>
            </div>
          </div>
          <div style={{background:"white",borderRadius:10,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <table style={{width:"100%",borderCollapse:"collapse",fontSize:13}}>
              <thead><tr style={{background:"#F9FAFB"}}>{["Sample","Events",gateLabel+" %","MFI ("+gateLabel+")","gMFI (all)","gMFI ("+gateLabel+")","log₂FC ("+basisLabel+")"].map(h=><th key={h} style={{padding:"10px 16px",textAlign:h==="Sample"?"left":"right",fontWeight:600,color:"#374151",borderBottom:"1px solid #E5E7EB",whiteSpace:"nowrap"}}>{h}</th>)}</tr></thead>
              <tbody>{rows.map((r,i)=>{
                const rowG=log2Basis==="all"?r.gAll:r.gGated;
                const isRef=r.s._idx===ref.s._idx;
                let l2;
                if(isRef)l2=<span style={{color:"#9CA3AF"}}>ref</span>;
                else if(rowG>0&&refG>0){const v=Math.log2(rowG/refG);l2=<span style={{fontWeight:700,color:v>0?"#15803D":v<0?"#B91C1C":"#374151"}}>{(v>0?"+":"")+v.toFixed(2)}</span>;}
                else l2=<span style={{color:"#9CA3AF"}}>—</span>;
                return <tr key={i} style={{borderBottom:i<rows.length-1?"1px solid #F3F4F6":"none",background:isRef?"#F9FAFB":"white"}}>
                  <td style={{padding:"9px 16px"}}><div style={{display:"flex",alignItems:"center",gap:8}}><div style={{width:10,height:10,borderRadius:3,background:colors[r.s._idx]||PALETTE[r.s._idx%PALETTE.length],flexShrink:0}}/><span style={{color:"#111827"}}>{r.s.name}</span></div></td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.n.toLocaleString()}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",fontWeight:700,color:gateColor}}>{r.pct}%</td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.mfiNum!=null?Math.round(r.mfiNum).toLocaleString():"—"}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.gAll!=null?Math.round(r.gAll).toLocaleString():"—"}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.gGated!=null?Math.round(r.gGated).toLocaleString():"—"}</td>
                  <td style={{padding:"9px 16px",textAlign:"right"}}>{l2}</td>
                </tr>;
              })}</tbody>
            </table>
          </div>
        </div>
        );
      })()}
    </>
  );
}

// ─── Quadrant Mode ───────────────────────────────────
function QuadrantMode({samples,allHeaders,colors,updateSampleName,updateColor,removeSample,swapSamples,moveSample,uploadProps,gateColor,showGate,showPct,initState,onState}){
  const S=initState||{};
  const[xCh,setXCh]=useState(S.xCh||"");
  const[yCh,setYCh]=useState(S.yCh||"");
  const[gateX,setGateX]=useState(S.gateX??500);
  const[gateY,setGateY]=useState(S.gateY??500);
  const[gateMode,setGateMode]=useState(S.gateMode||"quad"); // quad | vert | horiz
  const[dotMode,setDotMode]=useState(S.dotMode||"color"); // color | density
  const[xRangeIn,setXRangeIn]=useState(S.xRangeIn||{min:"",max:""});
  const[yRangeIn,setYRangeIn]=useState(S.yRangeIn||{min:"",max:""});
  const[xScale,setXScale]=useState(S.xScale||"log"); // log | biexp
  const[yScale,setYScale]=useState(S.yScale||"log");
  const[tickFont,setTickFont]=useState(S.tickFont??(S.axisFont||1));
  const[labelFont,setLabelFont]=useState(S.labelFont??(S.axisFont||1));
  const[dotSize,setDotSize]=useState(S.dotSize||2);
  const[showStats,setShowStats]=useState(S.showStats!==false);
  const[fcMode,setFcMode]=useState(S.fcMode||"gate"); // gate (− vs +) | sample (vs reference)
  const[fcRefIdx,setFcRefIdx]=useState(S.fcRefIdx!==undefined?S.fcRefIdx:null); // reference sample._idx
  const[fcBasis,setFcBasis]=useState(S.fcBasis||"all"); // all | pos — gMFI used for across-sample FC
  const[hidden,setHidden]=useState(S.hidden||[]); // sample._idx values hidden from the panel
  const toggleHide=idx=>setHidden(h=>h.includes(idx)?h.filter(i=>i!==idx):[...h,idx]);
  const[panelCols,setPanelCols]=useState(S.panelCols||"auto"); // "auto" | number of columns
  const[qDrag,setQDrag]=useState(null);   // sample._idx being dragged
  const[qOver,setQOver]=useState(null);   // sample._idx being dragged over
  const[qMenu,setQMenu]=useState(null);   // sample._idx whose ⋮ menu is open
  const qMenuRef=useRef(null);
  useEffect(()=>{if(qMenu===null)return;const h=e=>{if(qMenuRef.current&&!qMenuRef.current.contains(e.target))setQMenu(null);};document.addEventListener("mousedown",h);return()=>document.removeEventListener("mousedown",h);},[qMenu]);

  useEffect(()=>{if(onState)onState({xCh,yCh,gateX,gateY,gateMode,dotMode,xRangeIn,yRangeIn,xScale,yScale,tickFont,labelFont,dotSize,showStats,fcMode,fcRefIdx,fcBasis,hidden,panelCols});});

  useEffect(()=>{
    if(!allHeaders.length)return;
    const[dx,dy]=detectChannels(allHeaders);
    if(!xCh||!allHeaders.includes(xCh))setXCh(dx);
    if(!yCh||!allHeaders.includes(yCh))setYCh(dy);
  },[allHeaders,xCh,yCh]);

  const validSamples=useMemo(()=>
    samples.map((s,i)=>({...s,_idx:i})).filter(s=>
      s.columns[xCh]&&s.columns[xCh].length>0&&s.columns[yCh]&&s.columns[yCh].length>0
    ),[samples,xCh,yCh]);

  const xRangeAuto=useMemo(()=>{
    const pooled=[];
    for(const s of validSamples){const v=s.columns[xCh];if(v)for(let i=0;i<v.length;i++)pooled.push(v[i]);}
    if(!pooled.length)return{lo:0,hi:5,dMin:1,dMax:100000};
    return xScale==="log"?logRange(pooled):scatterRange(pooled);
  },[validSamples,xCh,xScale]);
  const yRangeAuto=useMemo(()=>{
    const pooled=[];
    for(const s of validSamples){const v=s.columns[yCh];if(v)for(let i=0;i<v.length;i++)pooled.push(v[i]);}
    if(!pooled.length)return{lo:0,hi:5,dMin:1,dMax:100000};
    return yScale==="log"?logRange(pooled):scatterRange(pooled);
  },[validSamples,yCh,yScale]);
  const xRange=xScale==="log"?manualLog(xRangeIn,xRangeAuto):manualBiexp(xRangeIn,xRangeAuto);
  const yRange=yScale==="log"?manualLog(yRangeIn,yRangeAuto):manualBiexp(yRangeIn,yRangeAuto);

  const shownSamples=useMemo(()=>validSamples.filter(s=>!hidden.includes(s._idx)),[validSamples,hidden]);
  const hiddenSamples=useMemo(()=>validSamples.filter(s=>hidden.includes(s._idx)),[validSamples,hidden]);
  const autoCols=shownSamples.length<=1?1:shownSamples.length<=4?2:3;
  const gridCols=panelCols==="auto"?autoCols:Math.max(1,Math.min(panelCols,Math.max(1,shownSamples.length)));
  const onGateChange=useCallback((x,y)=>{setGateX(x);setGateY(y);},[]);
  const fluoroHeaders=useMemo(()=>allHeaders.filter(h=>!/Time/i.test(h)),[allHeaders]);
  const stepBtn={padding:"3px 8px",border:"none",background:"white",fontSize:11,fontWeight:700,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"};

  // Composite the visible sample plots into ONE image so the whole panel exports as a single figure.
  const exportPanel=useCallback(()=>{
    if(!shownSamples.length)return;
    const {W,H}=DPG;const cols=gridCols,rows=Math.ceil(shownSamples.length/cols),scale=3;
    const c=document.createElement("canvas");c.width=W*cols*scale;c.height=H*rows*scale;
    const ctx=c.getContext("2d");ctx.scale(scale,scale);ctx.fillStyle="white";ctx.fillRect(0,0,W*cols,H*rows);
    shownSamples.forEach((s,i)=>{
      const cx=i%cols,cy=(i/cols)|0;
      paintDotPlot(ctx,cx*W,cy*H,{xVals:s.columns[xCh],yVals:s.columns[yCh],name:s.name,xCh,yCh,xRange,yRange,gateX,gateY,gateMode,dotMode,
        color:colors[s._idx]||PALETTE[s._idx%PALETTE.length],dotSize,tickFont,labelFont,showGate,showPct,gateColor,xScale,yScale,
        showYAxis:cx===0,showXAxis:i+cols>=shownSamples.length});
    });
    downloadDataURL(c.toDataURL("image/png"),"dot_plot_panel.png");
  },[shownSamples,gridCols,xCh,yCh,xRange,yRange,gateX,gateY,gateMode,dotMode,colors,dotSize,tickFont,labelFont,showGate,showPct,gateColor,xScale,yScale]);

  return(
    <>
      {/* Upload + Controls */}
      <div style={{maxWidth:1200,margin:"0 auto 20px",display:"flex",gap:14,flexWrap:"wrap",alignItems:"stretch"}}>
        <UploadBox {...uploadProps} samples={samples} accent="#0A9396" hint={<>One CSV per sample, <b>pre-gated to your population</b> (e.g. live singlets).<br/>Include both channels you'll plot on X &amp; Y (e.g. viability vs. signal).</>}/>
        {samples.length>0&&(
          <div style={{flex:"0 0 360px",background:"white",borderRadius:12,border:"1px solid #E5E7EB",padding:"14px 18px",display:"flex",flexDirection:"column",gap:12}}>
            <div style={{display:"flex",gap:12}}>
              <div style={{flex:1}}>
                <label style={labelStyle}>X-axis</label>
                <select value={xCh} onChange={e=>setXCh(e.target.value)} style={{...inputStyle,width:"100%",background:"white"}}>
                  {fluoroHeaders.map(h=><option key={h} value={h}>{h}</option>)}
                </select>
              </div>
              <div style={{flex:1}}>
                <label style={labelStyle}>Y-axis</label>
                <select value={yCh} onChange={e=>setYCh(e.target.value)} style={{...inputStyle,width:"100%",background:"white"}}>
                  {fluoroHeaders.map(h=><option key={h} value={h}>{h}</option>)}
                </select>
              </div>
            </div>
            <div style={{display:"flex",gap:12}}>
              {[["X scale",xScale,setXScale],["Y scale",yScale,setYScale]].map(([lbl,val,set])=>(
                <div key={lbl} style={{flex:1}}>
                  <label style={labelStyle}>{lbl}</label>
                  <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden",marginTop:3}}>
                    {[["log","Log"],["biexp","Biexp"]].map(([m,t],i)=>(
                      <button key={m} onClick={()=>set(m)} title={m==="log"?"Log scale (10ⁿ decades)":"Biexponential — spreads near-zero/negative events"} style={{flex:1,padding:"5px 4px",border:"none",borderLeft:i?"1px solid #E5E7EB":"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:val===m?"#EFF6FF":"white",color:val===m?"#3B82F6":"#9CA3AF"}}>{t}</button>
                    ))}
                  </div>
                </div>
              ))}
            </div>
            <div style={{display:"flex",gap:12}}>
              <div style={{flex:1}}>
                <label style={labelStyle}>Gate mode</label>
                <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden",marginTop:3}}>
                  {[["quad","Quad"],["vert","Vert"],["horiz","Horiz"]].map(([m,lbl],i)=>(
                    <button key={m} onClick={()=>setGateMode(m)} style={{flex:1,padding:"5px 4px",border:"none",borderLeft:i?"1px solid #E5E7EB":"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:gateMode===m?"#EFF6FF":"white",color:gateMode===m?"#3B82F6":"#9CA3AF"}}>{lbl}</button>
                  ))}
                </div>
              </div>
              <div style={{flex:1}}>
                <label style={labelStyle}>Dots</label>
                <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden",marginTop:3}}>
                  {[["color","Color"],["density","Density"]].map(([m,lbl],i)=>(
                    <button key={m} onClick={()=>setDotMode(m)} style={{flex:1,padding:"5px 4px",border:"none",borderLeft:i?"1px solid #E5E7EB":"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:dotMode===m?"#EFF6FF":"white",color:dotMode===m?"#3B82F6":"#9CA3AF"}}>{lbl}</button>
                  ))}
                </div>
              </div>
            </div>
            <div style={{display:"flex",gap:12}}>
              <div style={{flex:1}}>
                <label style={{...labelStyle,color:gateMode==="horiz"?"#D1D5DB":"#9CA3AF"}}>X gate</label>
                <input type="number" value={gateX} disabled={gateMode==="horiz"} onChange={e=>setGateX(Number(e.target.value))} style={{...inputStyle,width:"100%",background:gateMode==="horiz"?"#F9FAFB":"white",color:gateMode==="horiz"?"#9CA3AF":"#111827"}}/>
              </div>
              <div style={{flex:1}}>
                <label style={{...labelStyle,color:gateMode==="vert"?"#D1D5DB":"#9CA3AF"}}>Y gate</label>
                <input type="number" value={gateY} disabled={gateMode==="vert"} onChange={e=>setGateY(Number(e.target.value))} style={{...inputStyle,width:"100%",background:gateMode==="vert"?"#F9FAFB":"white",color:gateMode==="vert"?"#9CA3AF":"#111827"}}/>
              </div>
            </div>
            <div>
              <label style={labelStyle}>Axis range (raw units · blank = auto)</label>
              <AxisBar xLabel={xCh} yLabel={yCh} xR={xRangeIn} setXR={setXRangeIn} yR={yRangeIn} setYR={setYRangeIn} autoX={xRangeAuto} autoY={yRangeAuto}/>
            </div>
            <div style={{display:"flex",gap:10,flexWrap:"wrap"}}>
              {[["Tick font",tickFont,setTickFont,"Font size of the axis tick numbers"],["Label font",labelFont,setLabelFont,"Font size of the X/Y axis titles"]].map(([lbl,val,set,tip])=>(
                <div key={lbl} style={{flex:"1 1 100px"}} title={tip}>
                  <label style={labelStyle}>{lbl}</label>
                  <div style={{display:"flex",alignItems:"center",border:"1px solid #E5E7EB",borderRadius:6,overflow:"hidden",marginTop:3}}>
                    <button onClick={()=>set(f=>Math.max(0.6,+(f-0.1).toFixed(2)))} style={stepBtn}>A−</button>
                    <span style={{flex:1,fontSize:10,fontWeight:600,color:"#6B7280",textAlign:"center"}}>{Math.round(val*100)}%</span>
                    <button onClick={()=>set(f=>Math.min(2.2,+(f+0.1).toFixed(2)))} style={{...stepBtn,fontSize:13,borderLeft:"1px solid #E5E7EB"}}>A+</button>
                  </div>
                </div>
              ))}
              <div style={{flex:"1 1 100px"}} title="Weight (pixel size) of each event dot">
                <label style={labelStyle}>Dot size</label>
                <div style={{display:"flex",alignItems:"center",border:"1px solid #E5E7EB",borderRadius:6,overflow:"hidden",marginTop:3}}>
                  <button onClick={()=>setDotSize(d=>Math.max(1,d-1))} style={stepBtn}>−</button>
                  <span style={{flex:1,fontSize:10,fontWeight:600,color:"#6B7280",textAlign:"center"}}>{dotSize}px</span>
                  <button onClick={()=>setDotSize(d=>Math.min(6,d+1))} style={{...stepBtn,borderLeft:"1px solid #E5E7EB"}}>+</button>
                </div>
              </div>
            </div>
            <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",gap:8}}>
              <label style={{...labelStyle,margin:0}} title="Show the geometric MFI (whole vs gate) & log₂ fold-change table below, and include it in PNG exports">gMFI · log₂FC table</label>
              <button onClick={()=>setShowStats(s=>!s)} style={{padding:"4px 12px",borderRadius:6,border:"1px solid "+(showStats?"#3B82F6":"#E5E7EB"),background:showStats?"#EFF6FF":"white",color:showStats?"#3B82F6":"#9CA3AF",fontSize:11,fontWeight:700,cursor:"pointer",fontFamily:"var(--ff)"}}>{showStats?"On":"Off"}</button>
            </div>
            <div style={{display:"flex",alignItems:"center",gap:8,marginTop:"auto"}}>
              <span style={{fontSize:12,color:"#6B7280"}}>{validSamples.length+" sample"+(validSamples.length!==1?"s":"")}</span>
              <span style={{fontSize:11,color:"#9CA3AF"}}>{gateMode==="vert"?"Drag the vertical line":gateMode==="horiz"?"Drag the horizontal line":"Drag the crosshair"}</span>
            </div>
          </div>
        )}
      </div>

      {/* Dot Plots Panel */}
      {validSamples.length>0&&(
        <div style={{maxWidth:1200,margin:"0 auto",background:"white",borderRadius:12,border:"1px solid #E5E7EB",boxShadow:"0 1px 4px rgba(0,0,0,0.04)",padding:"14px 16px 16px"}}>
          <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",marginBottom:10}}>
            <span style={{fontSize:12,fontWeight:700,color:"#374151"}}>{shownSamples.length+" of "+validSamples.length+" shown · "+(gateMode==="quad"?"quadrant":gateMode)+" gate"}</span>
            <div style={{display:"flex",alignItems:"center",gap:8,marginLeft:"auto"}}>
              <span style={{fontSize:11,fontWeight:600,color:"#9CA3AF"}} title="Columns in the panel grid & export">Layout</span>
              <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
                {["auto",1,2,3,4,5].map((c,i)=>{const on=panelCols===c;return(
                  <button key={c} onClick={()=>setPanelCols(c)} title={c==="auto"?"Auto layout":c+(c===1?" column":" columns")}
                    style={{padding:"4px 9px",border:"none",borderLeft:i?"1px solid #E5E7EB":"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:on?"#EFF6FF":"white",color:on?"#3B82F6":"#9CA3AF"}}>{c==="auto"?"Auto":c}</button>
                );})}
              </div>
              <button onClick={exportPanel} disabled={!shownSamples.length} title="Export the whole panel as one PNG"
              style={{display:"flex",alignItems:"center",gap:6,padding:"6px 12px",borderRadius:7,border:"1px solid #E5E7EB",background:"white",color:shownSamples.length?"#374151":"#CBD5E1",fontSize:12,fontWeight:600,cursor:shownSamples.length?"pointer":"default",fontFamily:"var(--ff)"}}
              onMouseEnter={e=>{if(shownSamples.length)e.currentTarget.style.background="#F9FAFB";}} onMouseLeave={e=>{e.currentTarget.style.background="white";}}>
              <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
              Export panel
            </button>
            </div>
          </div>
          {hiddenSamples.length>0&&<div style={{display:"flex",alignItems:"center",gap:6,flexWrap:"wrap",marginBottom:10}}>
            <span style={{fontSize:11,fontWeight:600,color:"#9CA3AF"}}>Hidden:</span>
            {hiddenSamples.map(s=>(
              <button key={s._idx} onClick={()=>toggleHide(s._idx)} title="Show in panel"
                style={{display:"flex",alignItems:"center",gap:5,padding:"3px 9px",borderRadius:999,border:"1px solid #E5E7EB",background:"white",color:"#6B7280",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)"}}
                onMouseEnter={e=>{e.currentTarget.style.background="#F9FAFB";}} onMouseLeave={e=>{e.currentTarget.style.background="white";}}>
                <span style={{width:8,height:8,borderRadius:2,background:colors[s._idx]||PALETTE[s._idx%PALETTE.length],flexShrink:0}}/>
                {s.name}
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"/><circle cx="12" cy="12" r="3"/></svg>
              </button>
            ))}
            <button onClick={()=>setHidden([])} style={{padding:"3px 9px",borderRadius:999,border:"1px solid #E5E7EB",background:"white",color:"#6B7280",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)"}}>Show all</button>
          </div>}
          {shownSamples.length>0
            ?<div style={{display:"grid",gridTemplateColumns:"repeat("+gridCols+", 1fr)",gap:8}}>
            {shownSamples.map((s,pos)=>{
              const isOver=qOver===s._idx&&qDrag!==null&&qDrag!==s._idx;
              return(
              <div key={s._idx}
                onDragOver={moveSample&&shownSamples.length>1?e=>{e.preventDefault();setQOver(s._idx);}:undefined}
                onDragLeave={moveSample?()=>setQOver(o=>o===s._idx?null:o):undefined}
                onDrop={moveSample?e=>{e.preventDefault();if(qDrag!==null&&qDrag!==s._idx)moveSample(qDrag,s._idx);setQDrag(null);setQOver(null);}:undefined}
                style={{position:"relative",borderRadius:12,outline:isOver?"2px solid #3B82F6":"none",outlineOffset:2,opacity:qDrag===s._idx?0.5:1}}>
                {moveSample&&shownSamples.length>1&&<span draggable title="Drag to reorder"
                  onDragStart={e=>{setQDrag(s._idx);e.dataTransfer.effectAllowed="move";try{e.dataTransfer.setData("text/plain",String(s._idx));}catch(_){}}}
                  onDragEnd={()=>{setQDrag(null);setQOver(null);}}
                  style={{position:"absolute",top:4,left:4,zIndex:10,display:"flex",alignItems:"center",justifyContent:"center",width:26,height:26,borderRadius:6,background:"rgba(255,255,255,0.9)",border:"1px solid #E5E7EB",color:"#9CA3AF",cursor:"grab"}}><GripIcon/></span>}
                <div style={{position:"absolute",top:4,right:4,zIndex:10,display:"flex",gap:3}}>
                  <button onClick={()=>toggleHide(s._idx)} title="Hide from panel"
                    style={{width:26,height:26,borderRadius:6,border:"1px solid #E5E7EB",background:"rgba(255,255,255,0.9)",color:"#6B7280",cursor:"pointer",display:"flex",alignItems:"center",justifyContent:"center"}}
                    onMouseEnter={e=>{e.currentTarget.style.color="#3B82F6";}} onMouseLeave={e=>{e.currentTarget.style.color="#6B7280";}}><EyeIcon off={false}/></button>
                  <div style={{position:"relative"}} ref={qMenu===s._idx?qMenuRef:null}>
                    <button onClick={()=>setQMenu(qMenu===s._idx?null:s._idx)} title="More options"
                      style={{width:26,height:26,borderRadius:6,border:"1px solid #E5E7EB",background:qMenu===s._idx?"#F3F4F6":"rgba(255,255,255,0.9)",color:"#6B7280",cursor:"pointer",display:"flex",alignItems:"center",justifyContent:"center"}}><KebabIcon/></button>
                    {qMenu===s._idx&&<div style={{position:"absolute",top:"112%",right:0,zIndex:60,background:"white",border:"1px solid #E5E7EB",borderRadius:8,boxShadow:"0 8px 24px rgba(0,0,0,0.12)",padding:4,minWidth:130}}>
                      <button onClick={()=>{toggleHide(s._idx);setQMenu(null);}} style={{...menuItemStyle,color:"#374151"}} onMouseEnter={e=>e.currentTarget.style.background="#F3F4F6"} onMouseLeave={e=>e.currentTarget.style.background="none"}>Hide from panel</button>
                      <button onClick={()=>{setQMenu(null);if(window.confirm("Remove \""+s.name+"\" from the analysis?"))removeSample(s._idx);}} style={{...menuItemStyle,color:"#EF4444"}} onMouseEnter={e=>e.currentTarget.style.background="#FEF2F2"} onMouseLeave={e=>e.currentTarget.style.background="none"}>Delete sample</button>
                    </div>}
                  </div>
                </div>
                <DotPlot xVals={s.columns[xCh]} yVals={s.columns[yCh]} name={s.name} xCh={xCh} yCh={yCh}
                  xRange={xRange} yRange={yRange} gateX={gateX} gateY={gateY} onGateChange={onGateChange}
                  color={colors[s._idx]||PALETTE[s._idx%PALETTE.length]} onNameChange={n=>updateSampleName(s._idx,n)} gateColor={gateColor} showGate={showGate} showPct={showPct} gateMode={gateMode} dotMode={dotMode} xScale={xScale} yScale={yScale} tickFont={tickFont} labelFont={labelFont} dotSize={dotSize} showStats={showStats} bare/>
              </div>
              );
            })}
          </div>
            :<div style={{textAlign:"center",color:"#9CA3AF",fontSize:13,padding:"40px 0"}}>All plots hidden — use the chips above to show them.</div>}
        </div>
      )}

      {/* Summary table */}
      {validSamples.length>1&&(
        <div style={{maxWidth:1200,margin:"20px auto 0"}}>
          <div style={{background:"white",borderRadius:10,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <table style={{width:"100%",borderCollapse:"collapse",fontSize:13}}>
              <thead><tr style={{background:"#F9FAFB"}}>{(gateMode==="quad"?["Sample","Events","Q1 (−X/+Y)","Q2 (+X/+Y)","Q3 (−X/−Y)","Q4 (+X/−Y)"]:gateMode==="vert"?["Sample","Events",xCh+" − %",xCh+" + %"]:["Sample","Events",yCh+" + %",yCh+" − %"]).map(h=><th key={h} style={{padding:"10px 14px",textAlign:h==="Sample"?"left":"right",fontWeight:600,color:"#374151",borderBottom:"1px solid #E5E7EB",fontSize:12,whiteSpace:"nowrap"}}>{h}</th>)}</tr></thead>
              <tbody>{validSamples.map((s,i)=>{
                const xv=s.columns[xCh],yv=s.columns[yCh];
                const n=Math.min(xv.length,yv.length);
                let q1=0,q2=0,q3=0,q4=0,xp=0,yp=0;
                for(let k=0;k<n;k++){const xpos=xv[k]>=gateX,ypos=yv[k]>=gateY;if(xpos)xp++;if(ypos)yp++;if(!xpos&&ypos)q1++;else if(xpos&&ypos)q2++;else if(!xpos&&!ypos)q3++;else q4++;}
                const pc=v=>n?((v/n)*100).toFixed(1):"0.0";
                return(
                  <tr key={i} style={{borderBottom:i<validSamples.length-1?"1px solid #F3F4F6":"none"}}>
                    <td style={{padding:"9px 14px"}}><div style={{display:"flex",alignItems:"center",gap:8}}><div style={{width:10,height:10,borderRadius:3,background:colors[s._idx]||PALETTE[s._idx%PALETTE.length],flexShrink:0}}/><span style={{color:"#111827"}}>{s.name}</span></div></td>
                    <td style={{padding:"9px 14px",textAlign:"right",color:"#6B7280"}}>{n.toLocaleString()}</td>
                    {gateMode==="quad"?<>
                      <td style={{padding:"9px 14px",textAlign:"right",fontWeight:600,color:"#374151"}}>{pc(q1)}%</td>
                      <td style={{padding:"9px 14px",textAlign:"right",fontWeight:700,color:gateColor}}>{pc(q2)}%</td>
                      <td style={{padding:"9px 14px",textAlign:"right",color:"#6B7280"}}>{pc(q3)}%</td>
                      <td style={{padding:"9px 14px",textAlign:"right",fontWeight:600,color:"#374151"}}>{pc(q4)}%</td>
                    </>:gateMode==="vert"?<>
                      <td style={{padding:"9px 14px",textAlign:"right",color:"#6B7280"}}>{pc(n-xp)}%</td>
                      <td style={{padding:"9px 14px",textAlign:"right",fontWeight:700,color:gateColor}}>{pc(xp)}%</td>
                    </>:<>
                      <td style={{padding:"9px 14px",textAlign:"right",fontWeight:700,color:gateColor}}>{pc(yp)}%</td>
                      <td style={{padding:"9px 14px",textAlign:"right",color:"#6B7280"}}>{pc(n-yp)}%</td>
                    </>}
                  </tr>
                );
              })}</tbody>
            </table>
          </div>
        </div>
      )}

      {/* gMFI + log₂ fold-change table */}
      {showStats&&validSamples.length>0&&(()=>{
        // gMFI is a fluorescence readout of the channel the gate sits on: X for a vertical
        // gate, Y for a horizontal gate, and both for a quad gate — each split by its own gate.
        const readouts=gateMode==="vert"?[{ch:xCh,side:"x"}]:gateMode==="horiz"?[{ch:yCh,side:"y"}]:[{ch:xCh,side:"x"},{ch:yCh,side:"y"}];
        const gm=(vals,keep)=>{let s=0,c=0;for(let i=0;i<vals.length;i++){if((!keep||keep(i))&&vals[i]>0){s+=Math.log(vals[i]);c++;}}return c>0?Math.exp(s/c):null;};
        const fmt=v=>v==null?"—":Math.round(v).toLocaleString();
        const fcT=v=>v==null?"—":(v>0?"+":"")+v.toFixed(2);
        const l2=(a,b)=>(a>0&&b>0)?Math.log2(a/b):null;
        const posLbl=side=>side==="x"?xCh+"⁺":yCh+"⁺";
        // Per-sample, per-readout geometric MFIs.
        const data=validSamples.map(s=>{
          const xv=s.columns[xCh],yv=s.columns[yCh];
          return{s,per:readouts.map(r=>{
            const rv=r.ch===xCh?xv:yv;
            const keepPos=r.side==="x"?(k=>xv[k]>=gateX):(k=>yv[k]>=gateY);
            const keepNeg=r.side==="x"?(k=>xv[k]<gateX):(k=>yv[k]<gateY);
            return{all:gm(rv),pos:gm(rv,keepPos),neg:gm(rv,keepNeg)};
          })};
        });
        const ref=data.find(d=>d.s._idx===fcRefIdx)||data[0];
        const sampleMode=fcMode==="sample";
        const basisLbl=fcBasis==="all"?"all":"gate⁺";
        const fcHead=r=>sampleMode?"log₂FC "+r.ch+" vs "+ref.s.name:"log₂FC "+r.ch+" (+/−)";
        const head=["Sample"];
        readouts.forEach(r=>head.push(r.ch+" gMFI (all)",r.ch+" gMFI ("+posLbl(r.side)+")",fcHead(r)));
        return(
          <div style={{maxWidth:1200,margin:"20px auto 0"}}>
            <div style={{display:"flex",alignItems:"center",gap:10,flexWrap:"wrap",marginBottom:8,padding:"0 2px"}}>
              <span style={{fontSize:12,fontWeight:700,color:"#374151"}}>Geometric MFI &amp; log₂ fold-change</span>
              <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
                <button onClick={()=>setFcMode("gate")} title="log₂ ratio of the gate-positive gMFI to the gate-negative gMFI, within each sample" style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:!sampleMode?"#EFF6FF":"white",color:!sampleMode?"#3B82F6":"#9CA3AF"}}>Across gate (+/−)</button>
                <button onClick={()=>setFcMode("sample")} title="log₂ ratio of each sample's gMFI to a chosen reference sample" style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:sampleMode?"#EFF6FF":"white",color:sampleMode?"#3B82F6":"#9CA3AF"}}>Across samples</button>
              </div>
              {sampleMode&&<>
                <span style={{fontSize:11,color:"#9CA3AF"}}>vs</span>
                <select value={ref.s._idx} onChange={e=>setFcRefIdx(Number(e.target.value))} style={{...inputStyle,marginTop:0,padding:"4px 8px",background:"white",fontSize:12}}>
                  {data.map(d=><option key={d.s._idx} value={d.s._idx}>{d.s.name}</option>)}
                </select>
                <span style={{fontSize:11,color:"#9CA3AF"}}>using</span>
                <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
                  <button onClick={()=>setFcBasis("all")} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:fcBasis==="all"?"#EFF6FF":"white",color:fcBasis==="all"?"#3B82F6":"#9CA3AF"}}>gMFI (all)</button>
                  <button onClick={()=>setFcBasis("pos")} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:fcBasis==="pos"?"#EFF6FF":"white",color:fcBasis==="pos"?"#3B82F6":"#9CA3AF"}}>gMFI (gate⁺)</button>
                </div>
              </>}
              <span style={{fontSize:11,color:"#9CA3AF"}}>{sampleMode?"log₂FC = each sample's "+basisLbl+" gMFI ÷ "+ref.s.name+"'s":"log₂FC = gate-positive gMFI ÷ gate-negative gMFI"}</span>
            </div>
            <div style={{background:"white",borderRadius:10,border:"1px solid #E5E7EB",overflow:"auto"}}>
              <table style={{width:"100%",borderCollapse:"collapse",fontSize:13}}>
                <thead><tr style={{background:"#F9FAFB"}}>{head.map((h,j)=><th key={j} style={{padding:"10px 14px",textAlign:j===0?"left":"right",fontWeight:600,color:"#374151",borderBottom:"1px solid #E5E7EB",fontSize:12,whiteSpace:"nowrap"}}>{h}</th>)}</tr></thead>
                <tbody>{data.map((d,i)=>{
                  const isRef=sampleMode&&d.s._idx===ref.s._idx;
                  return(
                    <tr key={i} style={{borderBottom:i<data.length-1?"1px solid #F3F4F6":"none",background:isRef?"#F9FAFB":"white"}}>
                      <td style={{padding:"9px 14px"}}><div style={{display:"flex",alignItems:"center",gap:8}}><div style={{width:10,height:10,borderRadius:3,background:colors[d.s._idx]||PALETTE[d.s._idx%PALETTE.length],flexShrink:0}}/><span style={{color:"#111827"}}>{d.s.name}</span></div></td>
                      {d.per.map((p,ri)=>{
                        let fc,fcNode;
                        if(sampleMode){
                          if(isRef){fcNode=<span style={{color:"#9CA3AF"}}>ref</span>;}
                          else{const rp=ref.per[ri];fc=fcBasis==="all"?l2(p.all,rp.all):l2(p.pos,rp.pos);}
                        }else{fc=l2(p.pos,p.neg);}
                        if(fcNode===undefined)fcNode=<span style={{color:fc==null?"#9CA3AF":fc>0?"#15803D":fc<0?"#B91C1C":"#374151"}}>{fcT(fc)}</span>;
                        return[
                          <td key={ri+"a"} style={{padding:"9px 14px",textAlign:"right",color:"#6B7280"}}>{fmt(p.all)}</td>,
                          <td key={ri+"b"} style={{padding:"9px 14px",textAlign:"right",color:"#374151",fontWeight:600}}>{fmt(p.pos)}</td>,
                          <td key={ri+"c"} style={{padding:"9px 14px",textAlign:"right",fontWeight:700}}>{fcNode}</td>
                        ];
                      })}
                    </tr>
                  );
                })}</tbody>
              </table>
            </div>
          </div>
        );
      })()}
    </>
  );
}

// ─── Palette Picker ──────────────────────────────────
function PalettePicker({value,onChange}){
  return(
    <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
      <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em"}}>Palette</span>
      {Object.entries(PALETTES).map(([pname,cols])=>{
        const on=pname===value;
        return(
          <button key={pname} onClick={()=>onChange(pname)} title={pname}
            style={{display:"flex",alignItems:"center",padding:3,borderRadius:7,cursor:"pointer",
              border:"1.5px solid "+(on?"#111827":"#E5E7EB"),background:"white",transition:"all 0.12s"}}>
            <span style={{display:"flex",borderRadius:4,overflow:"hidden",boxShadow:"0 0 0 1px rgba(0,0,0,0.06)"}}>
              {cols.slice(0,6).map((c,i)=><span key={i} style={{width:11,height:16,background:c}}/>)}
            </span>
          </button>
        );
      })}
    </div>
  );
}

// ─── Gating (polygon) ────────────────────────────────
const AnalysisIcon=()=><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="M3 4h18l-7 8v6l-4 2v-8z"/></svg>;

// Gating math runs in LINEAR (raw) space so scatter plots look like FlowJo.
function pointInPoly(px,py,poly){
  let inside=false;
  for(let i=0,j=poly.length-1;i<poly.length;j=i++){
    const xi=poly[i].x,yi=poly[i].y,xj=poly[j].x,yj=poly[j].y;
    if(((yi>py)!==(yj>py))&&(px<(xj-xi)*(py-yi)/(yj-yi)+xi))inside=!inside;
  }
  return inside;
}
function gateIdx(sample,xCh,yCh,poly,parentIdx){
  const xv=sample.columns[xCh],yv=sample.columns[yCh];
  if(!xv||!yv||!poly||poly.length<3)return parentIdx?parentIdx.slice():[];
  const out=[];const src=parentIdx||null;const n=src?src.length:sample.tot;
  for(let k=0;k<n;k++){const e=src?src[k]:k;if(pointInPoly(xv[e],yv[e],poly))out.push(e);}
  return out;
}
function autoPoly(sample,xCh,yCh,parentIdx){
  const xv=sample.columns[xCh],yv=sample.columns[yCh];
  const src=parentIdx||null;const n=src?src.length:sample.tot;
  const xs=[],ys=[];const step=Math.max(1,Math.floor(n/5000));
  for(let k=0;k<n;k+=step){const e=src?src[k]:k;xs.push(xv[e]);ys.push(yv[e]);}
  xs.sort((a,b)=>a-b);ys.sort((a,b)=>a-b);
  const q=(a,p)=>a.length?a[Math.min(a.length-1,Math.floor(a.length*p))]:0;
  const cx=q(xs,0.5),cy=q(ys,0.5);
  const rx=Math.max((q(xs,0.9)-q(xs,0.1))*0.62,1);
  const ry=Math.max((q(ys,0.9)-q(ys,0.1))*0.62,1);
  const verts=[];const N=8;
  for(let i=0;i<N;i++){const a=Math.PI*2*i/N-Math.PI/2;verts.push({x:cx+Math.cos(a)*rx,y:cy+Math.sin(a)*ry});}
  return verts;
}
function linTicks(min,max){
  const span=max-min;if(!(span>0))return[min];
  const rough=span/5;const mag=10**Math.floor(Math.log10(rough));
  const step=([1,2,2.5,5,10].find(m=>rough<=m*mag)||10)*mag;
  const out=[];for(let v=Math.ceil(min/step)*step;v<=max+step*1e-6;v+=step)out.push(v);
  return out;
}
function fmtK(v){
  if(Math.abs(v)<1e-9)return "0";
  const a=Math.abs(v);
  if(a>=1e6)return (v/1e6).toFixed(1)+"M";
  if(a>=1e3)return Math.round(v/1e3)+"K";
  return String(Math.round(v));
}
// Jet-style density colormap LUT (blue → cyan → green → yellow → red)
const JET=(()=>{const cl=x=>Math.max(0,Math.min(1,x));const n=64;const lut=[];for(let i=0;i<n;i++){const t=i/(n-1);const r=cl(1.5-Math.abs(4*t-3)),g=cl(1.5-Math.abs(4*t-2)),b=cl(1.5-Math.abs(4*t-1));lut.push("rgb("+(r*255|0)+","+(g*255|0)+","+(b*255|0)+")");}return lut;})();
// Darker, saturated density map that reads well on white (no pastel mid-tones)
// Classic FACS "jet" rainbow. Blue is a thin floor for sparse events; green→yellow→red span the
// cluster densities so density differences read clearly (warm colors start early, not just at the peak).
const DENSITY=(()=>{const stops=["#0A1AC0","#1E5AF5","#12A8E8","#10D0C0","#22C82E","#8FD400","#E8D000","#FFB000","#FF7A00","#F84010","#D00000","#900000"].map(h=>[parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)]);const N=64;const lut=[];for(let i=0;i<N;i++){const t=i/(N-1)*(stops.length-1);const j=Math.min(stops.length-2,Math.floor(t));const f=t-j;const a=stops[j],b=stops[j+1];lut.push("rgb("+Math.round(a[0]+(b[0]-a[0])*f)+","+Math.round(a[1]+(b[1]-a[1])*f)+","+Math.round(a[2]+(b[2]-a[2])*f)+")");}return lut;})();

function GatePlot({title,sample,xCh,yCh,xRange,yRange,parentIdx,poly,onPoly,onCommit,onReset,gateColor,pct,count,parentCount}){
  const W=330,H=330,ML=52,MR=14,MT=26,MB=42,PW=W-ML-MR,PH=H-MT-MB;
  const canvasRef=useRef(null),svgRef=useRef(null),dragRef=useRef(null);
  const xmin=xRange.min,xmax=xRange.max,ymin=yRange.min,ymax=yRange.max;
  const xS=useCallback(v=>ML+((v-xmin)/(xmax-xmin))*PW,[xmin,xmax]);
  const yS=useCallback(v=>MT+PH-((v-ymin)/(ymax-ymin))*PH,[ymin,ymax]);
  const pxToVx=px=>xmin+((px-ML)/PW)*(xmax-xmin);
  const pyToVy=py=>ymin+((MT+PH-py)/PH)*(ymax-ymin);
  const xTicks=useMemo(()=>linTicks(xmin,xmax),[xmin,xmax]);
  const yTicks=useMemo(()=>linTicks(ymin,ymax),[ymin,ymax]);
  const gateName=title.split("·").pop().trim();

  useEffect(()=>{
    const canvas=canvasRef.current;if(!canvas)return;
    const ctx=canvas.getContext("2d");const dpr=window.devicePixelRatio||1;
    canvas.width=W*dpr;canvas.height=H*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,W,H);
    ctx.save();ctx.beginPath();ctx.rect(ML,MT,PW,PH);ctx.clip();
    const xv=sample.columns[xCh],yv=sample.columns[yCh];
    const src=parentIdx||null;const n=src?src.length:sample.tot;
    // 2-D density grid (pixel space) for heatmap coloring
    const BS=3,gc=Math.ceil(PW/BS),gr=Math.ceil(PH/BS);
    const grid=new Int32Array(gc*gr);
    const binOf=(a,b)=>{const gx=(a-ML)/BS|0,gy=(b-MT)/BS|0;return(gx<0||gy<0||gx>=gc||gy>=gr)?-1:gy*gc+gx;};
    for(let k=0;k<n;k++){const e=src?src[k]:k;const bi=binOf(xS(xv[e]),yS(yv[e]));if(bi>=0)grid[bi]++;}
    let maxD=1;for(let i=0;i<grid.length;i++)if(grid[i]>maxD)maxD=grid[i];
    const lm=Math.log(1+maxD);
    const maxPts=16000;const step=n>maxPts?Math.ceil(n/maxPts):1;
    ctx.globalAlpha=Math.max(0.5,Math.min(0.9,2600/Math.max(1,n)));
    for(let k=0;k<n;k+=step){
      const e=src?src[k]:k;const px=xS(xv[e]),py=yS(yv[e]);const bi=binOf(px,py);const d=bi>=0?grid[bi]:1;
      ctx.fillStyle=JET[Math.min(63,Math.floor(Math.log(1+d)/lm*63))];
      ctx.beginPath();ctx.arc(px,py,0.9,0,6.2832);ctx.fill();
    }
    ctx.restore();
  },[sample,xCh,yCh,xmin,xmax,ymin,ymax,parentIdx,xS,yS]);

  const svgPt=e=>{const r=svgRef.current.getBoundingClientRect();return{x:(e.clientX-r.left)/r.width*W,y:(e.clientY-r.top)/r.height*H};};
  const onDown=e=>{
    const{x,y}=svgPt(e);
    for(let i=0;i<poly.length;i++){if(Math.hypot(xS(poly[i].x)-x,yS(poly[i].y)-y)<9){dragRef.current={type:"vertex",i};e.preventDefault();return;}}
    if(pointInPoly(x,y,poly.map(p=>({x:xS(p.x),y:yS(p.y)})))){dragRef.current={type:"poly",sx:x,sy:y,orig:poly.map(p=>({x:p.x,y:p.y}))};e.preventDefault();}
  };
  const onMove=e=>{
    if(!dragRef.current)return;const{x,y}=svgPt(e);const d=dragRef.current;
    if(d.type==="vertex"){const np=poly.slice();np[d.i]={x:pxToVx(x),y:pyToVy(y)};onPoly(np);}
    else{const dvx=pxToVx(x)-pxToVx(d.sx),dvy=pyToVy(y)-pyToVy(d.sy);onPoly(d.orig.map(o=>({x:o.x+dvx,y:o.y+dvy})));}
  };
  const onUp=()=>{if(dragRef.current){dragRef.current=null;onCommit&&onCommit();}};
  const polyPath=poly.map((p,i)=>(i?"L":"M")+xS(p.x)+","+yS(p.y)).join(" ")+" Z";
  const cX=poly.reduce((a,p)=>a+xS(p.x),0)/poly.length;
  const cY=poly.reduce((a,p)=>a+yS(p.y),0)/poly.length;

  return(
    <div style={{background:"white",borderRadius:10,border:"1px solid #e2e5ea",padding:"8px 6px",boxShadow:"0 1px 4px rgba(0,0,0,0.04)"}}>
      <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",padding:"0 6px 4px"}}>
        <span style={{fontSize:12.5,fontWeight:700,color:"#111827"}}>{title}</span>
        <div style={{display:"flex",alignItems:"center",gap:8}}>
          <span style={{fontSize:12.5,fontWeight:800,color:gateColor}}>{pct}%</span>
          <button onClick={onReset} style={{...ghostBtn,padding:"3px 8px",fontSize:11}}>Reset</button>
        </div>
      </div>
      <div style={{position:"relative",width:W,maxWidth:"100%",aspectRatio:W+"/"+H,margin:"0 auto"}}>
        <canvas ref={canvasRef} style={{position:"absolute",top:0,left:0,width:"100%",height:"100%"}}/>
        <svg ref={svgRef} viewBox={"0 0 "+W+" "+H} style={{position:"absolute",top:0,left:0,width:"100%",height:"100%",userSelect:"none",touchAction:"none"}}
          onMouseDown={onDown} onMouseMove={onMove} onMouseUp={onUp} onMouseLeave={onUp}>
          <line x1={ML} x2={ML} y1={MT} y2={MT+PH} stroke="#4B5563" strokeWidth={1}/>
          <line x1={ML} x2={ML+PW} y1={MT+PH} y2={MT+PH} stroke="#111111" strokeWidth={1.2}/>
          {xTicks.map(t=>{const x=xS(t);if(x<ML-1||x>ML+PW+1)return null;return <g key={"x"+t}><line x1={x} x2={x} y1={MT+PH} y2={MT+PH+4} stroke="#111111"/><text x={x} y={MT+PH+15} textAnchor="middle" fontSize="8.5" fill="#4B5563" style={{fontFamily:"var(--ff)"}}>{fmtK(t)}</text></g>;})}
          {yTicks.map(t=>{const y=yS(t);if(y<MT-1||y>MT+PH+1)return null;return <g key={"y"+t}><line x1={ML-4} x2={ML} y1={y} y2={y} stroke="#4B5563"/><text x={ML-6} y={y+3} textAnchor="end" fontSize="8.5" fill="#4B5563" style={{fontFamily:"var(--ff)"}}>{fmtK(t)}</text></g>;})}
          <text x={ML+PW/2} y={H-4} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{xCh}</text>
          <text transform={"translate(12,"+(MT+PH/2)+") rotate(-90)"} textAnchor="middle" fontSize="10" fill="#374151" fontWeight="500" style={{fontFamily:"var(--ff)"}}>{yCh}</text>
          <path d={polyPath} fill="none" stroke="white" strokeWidth={3.4} strokeLinejoin="round" opacity={0.85}/>
          <path d={polyPath} fill="none" stroke={gateColor} strokeWidth={1.7} strokeLinejoin="round" style={{cursor:"move"}}/>
          {poly.map((p,i)=><circle key={i} cx={xS(p.x)} cy={yS(p.y)} r={4} fill="white" stroke={gateColor} strokeWidth={1.5} style={{cursor:"grab"}} onDoubleClick={ev=>{ev.stopPropagation();if(poly.length>3){onPoly(poly.filter((_,j)=>j!==i));onCommit&&onCommit();}}}/>)}
          <text x={cX} y={cY} textAnchor="middle" fontSize="11" fontWeight="700" fill="#111827" stroke="white" strokeWidth={3} paintOrder="stroke" style={{fontFamily:"var(--ff)",pointerEvents:"none"}}>{gateName}</text>
          <text x={cX} y={cY+13} textAnchor="middle" fontSize="10.5" fontWeight="700" fill={gateColor} stroke="white" strokeWidth={3} paintOrder="stroke" style={{fontFamily:"var(--ff)",pointerEvents:"none"}}>{pct}%</text>
          <text x={ML+PW-2} y={MT+11} textAnchor="end" fontSize="8" fill="#6B7280" style={{fontFamily:"var(--ff)"}}>{count.toLocaleString()+" / "+parentCount.toLocaleString()}</text>
        </svg>
      </div>
    </div>
  );
}

function AxisBar({xLabel,yLabel,xR,setXR,yR,setYR,autoX,autoY}){
  const box=(val,on,ph)=><input value={val} onChange={e=>on(e.target.value)} placeholder={ph} spellCheck={false}
    style={{width:56,padding:"3px 5px",borderRadius:5,border:"1px solid #D1D5DB",fontSize:10.5,fontFamily:"var(--ff)",color:"#111827",boxSizing:"border-box"}}/>;
  return <div style={{display:"flex",alignItems:"center",gap:5,flexWrap:"wrap",padding:"6px 8px 2px",fontSize:10.5,color:"#6B7280"}}>
    <span style={{fontWeight:700}}>{xLabel}</span>{box(xR.min,v=>setXR({...xR,min:v}),autoX?String(Math.round(autoX.dMin)):"min")}{box(xR.max,v=>setXR({...xR,max:v}),autoX?String(Math.round(autoX.dMax)):"max")}
    <span style={{fontWeight:700,marginLeft:6}}>{yLabel}</span>{box(yR.min,v=>setYR({...yR,min:v}),autoY?String(Math.round(autoY.dMin)):"min")}{box(yR.max,v=>setYR({...yR,max:v}),autoY?String(Math.round(autoY.dMax)):"max")}
    <button onClick={()=>{setXR({min:"",max:""});setYR({min:"",max:""});}} style={{...ghostBtn,padding:"2px 7px",fontSize:10.5,marginLeft:"auto"}}>Auto</button>
  </div>;
}

// ─── Analysis Mode (FCS → gating → plots) ────────────
function AnalysisMode({gateColor,showGate,showPct,onToggleGate,onTogglePct,paletteName}){
  const[fsamples,setFsamples]=useState([]);
  const[dispIdx,setDispIdx]=useState(0);
  const[fscA,setFscA]=useState("");
  const[sscA,setSscA]=useState("");
  const[fscH,setFscH]=useState("");
  const[cellsPoly,setCellsPoly]=useState(null);
  const[singletsPoly,setSingletsPoly]=useState(null);
  const[appCells,setAppCells]=useState(null);
  const[appSing,setAppSing]=useState(null);
  const[downCh,setDownCh]=useState("");
  const[downGate,setDownGate]=useState(1000);
  const[downView,setDownView]=useState("overlay");
  const[downNorm,setDownNorm]=useState(true);
  const[overlayPctPos,setOverlayPctPos]=useState("left");
  const[dragOver,setDragOver]=useState(false);
  const[busy,setBusy]=useState(false);
  const[cellsXR,setCellsXR]=useState({min:"",max:""});
  const[cellsYR,setCellsYR]=useState({min:"",max:""});
  const[singXR,setSingXR]=useState({min:"",max:""});
  const[singYR,setSingYR]=useState({min:"",max:""});
  const fileRef=useRef(null);

  const pal=PALETTES[paletteName]||PALETTES.Classic;
  const colors=useMemo(()=>fsamples.map((_,i)=>pal[i%pal.length]),[fsamples,pal]);

  const readBuf=file=>new Promise((res,rej)=>{const r=new FileReader();r.onload=()=>res(r.result);r.onerror=()=>rej(r.error);r.readAsArrayBuffer(file);});
  const handleFiles=async files=>{
    setBusy(true);const ns=[];
    for(const file of files){
      try{const buf=await readBuf(file);const p=parseFCS(buf);ns.push({name:file.name.replace(/\.fcs$/i,""),headers:p.headers,columns:p.columns,longMap:p.longMap,ranges:p.ranges,tot:p.tot});}
      catch(e){console.error("FCS parse failed:",file.name,e);}
    }
    setBusy(false);
    if(!ns.length)return;
    setFsamples(prev=>{
      if(prev.length===0){
        const sc=detectScatter(ns[0].headers);
        setFscA(sc.fscA||ns[0].headers[0]);setSscA(sc.sscA||ns[0].headers[1]||ns[0].headers[0]);setFscH(sc.fscH||ns[0].headers[0]);
        setDownCh(ns[0].headers.find(h=>!/FSC|SSC|Time|-H$|-W$|\.\d/i.test(h))||ns[0].headers[0]);
      }
      return[...prev,...ns];
    });
  };
  const onDrop=e=>{e.preventDefault();setDragOver(false);const f=[...e.dataTransfer.files].filter(f=>/\.fcs$/i.test(f.name));if(f.length)handleFiles(f);};
  const clearAll=()=>{setFsamples([]);setCellsPoly(null);setSingletsPoly(null);setAppCells(null);setAppSing(null);setDispIdx(0);};

  // initialise polygons once samples + channels ready
  useEffect(()=>{
    if(!fsamples.length||!fscA||!sscA||!fscH||cellsPoly)return;
    const s=fsamples[0];
    const cp=autoPoly(s,fscA,sscA,null);
    const cIdx=gateIdx(s,fscA,sscA,cp,null);
    const sp=autoPoly(s,fscA,fscH,cIdx);
    setCellsPoly(cp);setAppCells(cp);setSingletsPoly(sp);setAppSing(sp);
  },[fsamples,fscA,sscA,fscH,cellsPoly]);

  const disp=fsamples[dispIdx]||null;
  const dispGate=useMemo(()=>{
    if(!disp||!cellsPoly||!singletsPoly)return null;
    const cIdx=gateIdx(disp,fscA,sscA,cellsPoly,null);
    const cSet=new Uint8Array(disp.tot);for(const e of cIdx)cSet[e]=1;
    const sIdx=gateIdx(disp,fscA,fscH,singletsPoly,cIdx);
    const sSet=new Uint8Array(disp.tot);for(const e of sIdx)sSet[e]=1;
    return{cIdx,cSet,sIdx,sSet};
  },[disp,fscA,sscA,fscH,cellsPoly,singletsPoly]);

  // Full-scale LINEAR ranges (0 → instrument max) — consistent across both plots, shows every event
  const chRange=useCallback(ch=>{
    if(!disp)return null;
    let mx=disp.ranges&&disp.ranges[ch]?disp.ranges[ch]:0;
    if(!mx){const col=disp.columns[ch];if(col){const st=Math.max(1,col.length/5000|0);for(let i=0;i<col.length;i+=st)if(col[i]>mx)mx=col[i];mx*=1.05;}else mx=262144;}
    return{min:0,max:mx,dMin:0,dMax:mx};
  },[disp]);
  const cRX=useMemo(()=>chRange(fscA),[chRange,fscA]);
  const cRY=useMemo(()=>chRange(sscA),[chRange,sscA]);
  const sRX=useMemo(()=>chRange(fscA),[chRange,fscA]);
  const sRY=useMemo(()=>chRange(fscH),[chRange,fscH]);
  const mrange=(r,auto)=>{if(!auto)return auto;const mn=Number(r.min),mx=Number(r.max);if(!r.min.trim()||!r.max.trim()||!isFinite(mn)||!isFinite(mx)||mx<=mn)return auto;return{min:mn,max:mx,dMin:mn,dMax:mx};};
  const cxR=mrange(cellsXR,cRX),cyR=mrange(cellsYR,cRY),sxR=mrange(singXR,sRX),syR=mrange(singYR,sRY);

  const stats=useMemo(()=>{
    if(!appCells||!appSing||!fscA||!sscA||!fscH)return[];
    return fsamples.map(s=>{
      const cIdx=gateIdx(s,fscA,sscA,appCells,null);
      const sIdx=gateIdx(s,fscA,fscH,appSing,cIdx);
      return{name:s.name,tot:s.tot,cellsN:cIdx.length,singN:sIdx.length,singIdx:sIdx};
    });
  },[fsamples,fscA,sscA,fscH,appCells,appSing]);

  const gatedSamples=useMemo(()=>{
    if(!downCh)return[];
    return fsamples.map((s,i)=>{
      const idx=stats[i]?stats[i].singIdx:[];const col=s.columns[downCh];
      const vals=new Float32Array(idx.length);if(col)for(let k=0;k<idx.length;k++)vals[k]=col[idx[k]];
      return{name:s.name,columns:{[downCh]:vals}};
    });
  },[fsamples,stats,downCh]);

  const headers=disp?disp.headers:[];
  const fluoroHeaders=headers.filter(h=>!/Time/i.test(h));
  const chLabel=h=>disp&&disp.longMap[h]&&disp.longMap[h]!==h?h+" · "+disp.longMap[h]:h;
  const selStyle={...inputStyle,width:"100%",background:"white"};

  return(
    <>
      {/* Upload + scatter channel controls */}
      <div style={{maxWidth:1200,margin:"0 auto 16px",display:"flex",gap:14,flexWrap:"wrap",alignItems:"stretch"}}>
        <div onDragOver={e=>{e.preventDefault();setDragOver(true);}} onDragLeave={()=>setDragOver(false)} onDrop={onDrop}
          onClick={()=>fileRef.current&&fileRef.current.click()}
          style={{flex:"1 1 260px",minHeight:96,borderRadius:12,border:"2px dashed "+(dragOver?"#7C3AED":"#D1D5DB"),background:dragOver?"#F5F3FF":"white",display:"flex",flexDirection:"column",alignItems:"center",justifyContent:"center",cursor:"pointer",transition:"all 0.15s",padding:"14px 20px"}}>
          <svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke={dragOver?"#7C3AED":"#9CA3AF"} strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg>
          <span style={{fontSize:13,fontWeight:600,color:dragOver?"#7C3AED":"#374151",marginTop:5}}>{busy?"Reading…":(fsamples.length>0?"+ Add more .fcs files":"Drop raw .fcs files or click to browse")}</span>
          <span style={{fontSize:11,color:"#9CA3AF",marginTop:2}}>Ungated FCS from your cytometer — gating happens here</span>
          <input ref={fileRef} type="file" accept=".fcs" multiple onChange={e=>{handleFiles([...e.target.files]);e.target.value="";}} style={{display:"none"}}/>
        </div>
        {fsamples.length>0&&(
          <div style={{flex:"0 0 380px",background:"white",borderRadius:12,border:"1px solid #E5E7EB",padding:"14px 18px",display:"flex",flexDirection:"column",gap:10}}>
            <div>
              <label style={labelStyle}>Displayed sample</label>
              <select value={dispIdx} onChange={e=>setDispIdx(Number(e.target.value))} style={selStyle}>
                {fsamples.map((s,i)=><option key={i} value={i}>{s.name}</option>)}
              </select>
            </div>
            <div style={{display:"flex",gap:10}}>
              <div style={{flex:1}}><label style={labelStyle}>Cells X</label><select value={fscA} onChange={e=>setFscA(e.target.value)} style={selStyle}>{headers.map(h=><option key={h} value={h}>{h}</option>)}</select></div>
              <div style={{flex:1}}><label style={labelStyle}>Cells Y</label><select value={sscA} onChange={e=>setSscA(e.target.value)} style={selStyle}>{headers.map(h=><option key={h} value={h}>{h}</option>)}</select></div>
              <div style={{flex:1}}><label style={labelStyle}>Singlet Y</label><select value={fscH} onChange={e=>setFscH(e.target.value)} style={selStyle}>{headers.map(h=><option key={h} value={h}>{h}</option>)}</select></div>
            </div>
            <div style={{display:"flex",alignItems:"center",gap:8,marginTop:"auto"}}>
              <span style={{fontSize:12,color:"#6B7280"}}>{fsamples.length+" file"+(fsamples.length!==1?"s":"")}</span>
              <span style={{fontSize:11,color:"#9CA3AF"}}>Drag polygon vertices · double-click to remove</span>
              <button onClick={clearAll} style={{marginLeft:"auto",...ghostBtn}}>Clear all</button>
            </div>
          </div>
        )}
      </div>

      {/* Gating pipeline */}
      {disp&&cellsPoly&&singletsPoly&&dispGate&&cxR&&cyR&&sxR&&syR&&(
        <div style={{maxWidth:1200,margin:"0 auto 16px",display:"grid",gridTemplateColumns:"repeat(auto-fit,minmax(300px,1fr))",gap:16}}>
          <div>
            <GatePlot title="1 · Cells" sample={disp} xCh={fscA} yCh={sscA} xRange={cxR} yRange={cyR} parentIdx={null} gatedSet={dispGate.cSet}
              poly={cellsPoly} onPoly={setCellsPoly} onCommit={()=>setAppCells(cellsPoly)} onReset={()=>{const p=autoPoly(disp,fscA,sscA,null);setCellsPoly(p);setAppCells(p);}}
              color={colors[dispIdx]||pal[0]} gateColor={gateColor} pct={disp.tot?((dispGate.cIdx.length/disp.tot)*100).toFixed(1):"0.0"} count={dispGate.cIdx.length} parentCount={disp.tot}/>
            <AxisBar xLabel={fscA} yLabel={sscA} xR={cellsXR} setXR={setCellsXR} yR={cellsYR} setYR={setCellsYR} autoX={cRX} autoY={cRY}/>
          </div>
          <div>
            <GatePlot title="2 · Singlets" sample={disp} xCh={fscA} yCh={fscH} xRange={sxR} yRange={syR} parentIdx={dispGate.cIdx} gatedSet={dispGate.sSet}
              poly={singletsPoly} onPoly={setSingletsPoly} onCommit={()=>setAppSing(singletsPoly)} onReset={()=>{const p=autoPoly(disp,fscA,fscH,dispGate.cIdx);setSingletsPoly(p);setAppSing(p);}}
              color={colors[dispIdx]||pal[0]} gateColor={gateColor} pct={dispGate.cIdx.length?((dispGate.sIdx.length/dispGate.cIdx.length)*100).toFixed(1):"0.0"} count={dispGate.sIdx.length} parentCount={dispGate.cIdx.length}/>
            <AxisBar xLabel={fscA} yLabel={fscH} xR={singXR} setXR={setSingXR} yR={singYR} setYR={setSingYR} autoX={sRX} autoY={sRY}/>
          </div>
        </div>
      )}

      {/* Population table */}
      {stats.length>0&&(
        <div style={{maxWidth:1200,margin:"0 auto 18px"}}>
          <div style={{background:"white",borderRadius:10,border:"1px solid #E5E7EB",overflow:"hidden"}}>
            <table style={{width:"100%",borderCollapse:"collapse",fontSize:13}}>
              <thead><tr style={{background:"#F9FAFB"}}>{["Sample","Events","Cells","Cells %","Singlets","Singlets % (of cells)"].map(h=><th key={h} style={{padding:"10px 16px",textAlign:h==="Sample"?"left":"right",fontWeight:600,color:"#374151",borderBottom:"1px solid #E5E7EB",whiteSpace:"nowrap"}}>{h}</th>)}</tr></thead>
              <tbody>{stats.map((r,i)=>(
                <tr key={i} style={{borderBottom:i<stats.length-1?"1px solid #F3F4F6":"none"}}>
                  <td style={{padding:"9px 16px"}}><div style={{display:"flex",alignItems:"center",gap:8}}><div style={{width:10,height:10,borderRadius:3,background:colors[i]||pal[0]}}/><span style={{color:"#111827"}}>{r.name}</span></div></td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.tot.toLocaleString()}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.cellsN.toLocaleString()}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",fontWeight:600,color:"#374151"}}>{r.tot?((r.cellsN/r.tot)*100).toFixed(1):"0.0"}%</td>
                  <td style={{padding:"9px 16px",textAlign:"right",color:"#6B7280"}}>{r.singN.toLocaleString()}</td>
                  <td style={{padding:"9px 16px",textAlign:"right",fontWeight:700,color:gateColor}}>{r.cellsN?((r.singN/r.cellsN)*100).toFixed(1):"0.0"}%</td>
                </tr>
              ))}</tbody>
            </table>
          </div>
        </div>
      )}

      {/* Downstream plots of the singlet population */}
      {stats.length>0&&(
        <div style={{maxWidth:1200,margin:"0 auto"}}>
          <div style={{display:"flex",alignItems:"center",gap:10,flexWrap:"wrap",marginBottom:10}}>
            <span style={{fontSize:13,fontWeight:700,color:"#111827"}}>Singlet population —</span>
            <select value={downCh} onChange={e=>setDownCh(e.target.value)} style={{...inputStyle,marginTop:0,width:"auto",minWidth:160,background:"white"}}>
              {fluoroHeaders.map(h=><option key={h} value={h}>{chLabel(h)}</option>)}
            </select>
            <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden",marginLeft:4}}>
              <button onClick={()=>setDownView("overlay")} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:downView==="overlay"?"#EFF6FF":"white",color:downView==="overlay"?"#3B82F6":"#9CA3AF"}}>Overlay</button>
              <button onClick={()=>setDownView("grid")} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:downView==="grid"?"#EFF6FF":"white",color:downView==="grid"?"#3B82F6":"#9CA3AF"}}>Grid</button>
            </div>
            <span style={{fontSize:11,color:"#9CA3AF",marginLeft:4}}>Drag the gate to set positive %</span>
          </div>
          {downView==="overlay"
            ?<OverlayHistogram samples={gatedSamples} colors={colors} channel={downCh} xLabel={downCh} yLabel="Count" gateValue={downGate} onGateChange={setDownGate}
                xDomain={null} gateLabel="+" normalize={downNorm} gateColor={gateColor} showGate={showGate} showPct={showPct} onToggleGate={onToggleGate} onTogglePct={onTogglePct} pctPos={overlayPctPos} setPctPos={setOverlayPctPos}/>
            :<div style={{display:"grid",gridTemplateColumns:"repeat("+(gatedSamples.length===1?1:gatedSamples.length<=4?2:3)+", 1fr)",gap:16}}>
              {gatedSamples.map((s,i)=><Histogram key={i} values={s.columns[downCh]} name={s.name} color={colors[i]||pal[0]} xLabel={downCh} yLabel="Count"
                gateValue={downGate} onGateChange={setDownGate} gateLabel="+" gateColor={gateColor} showGate={showGate} showPct={showPct}/>)}
            </div>}
        </div>
      )}

      {fsamples.length===0&&<div style={{maxWidth:620,margin:"36px auto",textAlign:"center",color:"#9CA3AF",fontSize:13,lineHeight:1.7}}><strong style={{color:"#6B7280"}}>Full-stack analysis:</strong><br/>Drop raw <b>.fcs</b> files → adjust the Cells and Singlets polygon gates → plot any channel on the gated singlets.</div>}
    </>
  );
}

// ─── Shared button styles + Modal ────────────────────
const primaryBtn={padding:"7px 12px",borderRadius:7,border:"1px solid #2563EB",background:"#2563EB",color:"white",fontSize:12.5,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)"};
const ghostBtn={padding:"6px 11px",borderRadius:7,border:"1px solid #E5E7EB",background:"white",color:"#374151",fontSize:12,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)"};

function Modal({title,onClose,children,width=400}){
  return(
    <div onClick={onClose} style={{position:"fixed",inset:0,background:"rgba(17,24,39,0.45)",display:"flex",alignItems:"center",justifyContent:"center",zIndex:1000,padding:16}}>
      <div onClick={e=>e.stopPropagation()} style={{background:"white",borderRadius:14,border:"1px solid #E5E7EB",padding:20,width,maxWidth:"92vw",maxHeight:"86vh",overflow:"auto",boxShadow:"0 20px 60px rgba(0,0,0,0.25)",fontFamily:"var(--ff)"}}>
        <div style={{display:"flex",justifyContent:"space-between",alignItems:"center",marginBottom:14}}>
          <span style={{fontSize:16,fontWeight:700,color:"#111827"}}>{title}</span>
          <button onClick={onClose} style={{border:"none",background:"none",cursor:"pointer",fontSize:22,color:"#9CA3AF",lineHeight:1}}>×</button>
        </div>
        {children}
      </div>
    </div>
  );
}

function MergeChannelsModal({headers,onMerge,onClose}){
  const[sel,setSel]=useState([]);
  const[canon,setCanon]=useState("");
  const toggle=h=>setSel(s=>{const n=s.includes(h)?s.filter(x=>x!==h):[...s,h];return n;});
  useEffect(()=>{if(sel.length&&!sel.includes(canon))setCanon(sel[0]);if(!sel.length)setCanon("");},[sel]);
  return(
    <Modal title="Merge channels" onClose={onClose} width={430}>
      <p style={{fontSize:12.5,color:"#6B7280",margin:"0 0 12px",lineHeight:1.5}}>Pick channels that are really the same detector under different names (e.g. <b>mCherry</b> and <b>YL2-A</b>), then merge them into one so those samples plot together.</p>
      <div style={{maxHeight:250,overflowY:"auto",border:"1px solid #E5E7EB",borderRadius:8,padding:6,display:"flex",flexDirection:"column",gap:2}}>
        {headers.map(h=><label key={h} style={{display:"flex",alignItems:"center",gap:8,padding:"5px 7px",borderRadius:6,cursor:"pointer",fontSize:13,color:"#374151",background:sel.includes(h)?"#EFF6FF":"transparent"}}>
          <input type="checkbox" checked={sel.includes(h)} onChange={()=>toggle(h)}/>{h}
        </label>)}
      </div>
      {sel.length>=2
        ?<div style={{marginTop:12,display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
          <span style={{fontSize:12.5,fontWeight:600,color:"#374151"}}>Merge into:</span>
          <input value={canon} onChange={e=>setCanon(e.target.value)} style={{...inputStyle,marginTop:0,flex:1,minWidth:120}}/>
          <button onClick={()=>{onMerge(sel,canon.trim()||sel[0]);onClose();}} style={primaryBtn}>Merge</button>
        </div>
        :<div style={{marginTop:10,fontSize:11.5,color:"#9CA3AF"}}>Select at least two channels to merge. (Tip: name the result something clear like “mCherry / YL2-A”.)</div>}
    </Modal>
  );
}

// ─── Main App ────────────────────────────────────────
export default function FlowCytometryApp(){
  const[samples,setSamples]=useState([]);
  const[colors,setColors]=useState([]);
  const[dragOver,setDragOver]=useState(false);
  const[mode,setMode]=useState("histogram"); // "histogram" | "quadrant"
  const[paletteName,setPaletteName]=useState("Classic");
  const[gateColor,setGateColor]=useState("#2563EB");
  const[showGate,setShowGate]=useState(true);
  const[showPct,setShowPct]=useState(true);
  const[hidden,setHidden]=useState([]); // sample positions hidden from plots
  const[mergeOpen,setMergeOpen]=useState(false);
  const fileRef=useRef(null);
  // derived: union of all channel names across loaded samples (updates on merge)
  const allHeaders=useMemo(()=>{const set=new Set();for(const s of samples)for(const h of(s.headers||[]))set.add(h);return[...set];},[samples]);

  useEffect(()=>{try{const link=document.createElement("link");link.href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700;800&display=swap";link.rel="stylesheet";document.head.appendChild(link);}catch(e){}},[]);

  const readFileAsText=file=>new Promise((resolve,reject)=>{
    const reader=new FileReader();
    reader.onload=()=>resolve(reader.result);
    reader.onerror=()=>reject(reader.error);
    reader.readAsText(file);
  });
  const handleFiles=async files=>{
    const ns=[];
    for(let fi=0;fi<files.length;fi++){
      try{
        const file=files[fi];
        const text=await readFileAsText(file);
        const{headers,columns}=parseCSV(text);
        ns.push({name:file.name.replace(/\.(csv|tsv|txt)$/i,""),headers,columns});
      }catch(e){console.error("Failed to read file",e);}
    }
    const pal=PALETTES[paletteName]||PALETTES.Classic;
    setSamples(prev=>{const next=[...prev,...ns];setColors(pc=>{const nc=[...pc];for(let i=pc.length;i<next.length;i++)nc.push(pal[i%pal.length]);return nc;});return next;});
  };
  const onDrop=e=>{e.preventDefault();setDragOver(false);const f=[...e.dataTransfer.files].filter(f=>/\.(csv|tsv|txt)$/i.test(f.name));if(f.length)handleFiles(f);};
  const updateSampleName=(idx,n)=>setSamples(p=>p.map((s,i)=>i===idx?{...s,name:n}:s));
  const updateColor=(idx,c)=>setColors(p=>{const nc=[...p];nc[idx]=c;return nc;});
  const removeSample=idx=>{
    setSamples(p=>p.filter((_,i)=>i!==idx));
    setColors(p=>p.filter((_,i)=>i!==idx));
    setHidden(h=>h.filter(i=>i!==idx).map(i=>i>idx?i-1:i));
  };
  const toggleHidden=idx=>setHidden(h=>h.includes(idx)?h.filter(i=>i!==idx):[...h,idx]);
  const swapSamples=(i,j)=>{
    setSamples(p=>{if(i<0||j<0||i>=p.length||j>=p.length)return p;const a=[...p];[a[i],a[j]]=[a[j],a[i]];return a;});
    setColors(p=>{if(i<0||j<0||i>=p.length||j>=p.length)return p;const a=[...p];[a[i],a[j]]=[a[j],a[i]];return a;});
    setHidden(h=>h.map(x=>x===i?j:x===j?i:x));
  };
  // Move sample from index `from` to index `to` (drag-and-drop reorder), keeping colors + hidden in sync.
  const moveSample=(from,to)=>{
    if(from===to||from<0||to<0||from>=samples.length||to>=samples.length)return;
    const mv=arr=>{const a=[...arr];const[x]=a.splice(from,1);a.splice(to,0,x);return a;};
    setSamples(mv);setColors(mv);
    setHidden(h=>h.map(idx=>idx===from?to:from<to?(idx>from&&idx<=to?idx-1:idx):(idx>=to&&idx<from?idx+1:idx)));
  };
  const clearAll=()=>{setSamples([]);setColors([]);setHidden([]);};
  // Unify differently-named channels that are the same detector (e.g. "mCherry" and "YL2-A")
  const mergeChannels=(nameArr,canonical)=>{
    const names=nameArr.filter(Boolean);if(names.length<1||!canonical)return;
    setSamples(prev=>prev.map(s=>{
      const cols={...s.columns};
      for(const n of names){if(n===canonical)continue;if(cols[n]!==undefined){if(cols[canonical]===undefined)cols[canonical]=cols[n];delete cols[n];}}
      const seen=new Set();const nh=[];
      for(const h of(s.headers||[])){const m=names.includes(h)?canonical:h;if(!seen.has(m)){seen.add(m);nh.push(m);}}
      return{...s,columns:cols,headers:nh};
    }));
  };
  const applyPalette=name=>{setPaletteName(name);const cols=PALETTES[name]||PALETTES.Classic;setColors(samples.map((_,i)=>cols[i%cols.length]));};

  const uploadProps={dragOver,setDragOver,onDrop,fileRef,handleFiles};
  const accent=mode==="histogram"?"#3B82F6":"#0A9396";

  const modeBtn=(id,icon,label)=>{
    const on=mode===id;
    return <button onClick={()=>setMode(id)} style={{display:"flex",alignItems:"center",gap:6,padding:"7px 14px",border:"none",borderRadius:8,cursor:"pointer",fontFamily:"var(--ff)",fontSize:13,fontWeight:700,background:on?"white":"transparent",color:on?"#111827":"#6B7280",boxShadow:on?"0 1px 3px rgba(0,0,0,0.1)":"none",transition:"all 0.12s"}}>{icon} {label}</button>;
  };

  return(
    <div style={{"--ff":"'IBM Plex Sans', system-ui, sans-serif",fontFamily:"var(--ff)",minHeight:"100vh",background:"#F8F9FB",padding:"28px 24px"}}>

      {/* Header */}
      <div style={{maxWidth:1200,margin:"0 auto 20px"}}>
        <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",gap:16,flexWrap:"wrap"}}>
          <div>
            <div style={{display:"flex",alignItems:"center",gap:10,marginBottom:4}}>
              <div style={{width:32,height:32,borderRadius:8,background:mode==="histogram"?"linear-gradient(135deg,#3B82F6,#8B5CF6)":mode==="quadrant"?"linear-gradient(135deg,#0A9396,#CA391D)":"linear-gradient(135deg,#7C3AED,#DB2777)",display:"flex",alignItems:"center",justifyContent:"center",color:"white"}}>
                {mode==="histogram"?<HistIcon/>:mode==="quadrant"?<QuadIcon/>:<AnalysisIcon/>}
              </div>
              <h1 style={{fontSize:22,fontWeight:700,color:"#111827",margin:0}}>Flume</h1>
            </div>
            <p style={{color:"#6B7280",fontSize:13,margin:"4px 0 0 42px"}}>
              {mode==="histogram"
                ?"Click titles to rename · Click legend swatches to recolor · Drag the gate line · Overlay to superimpose samples"
                :mode==="quadrant"
                ?"Upload gated singlet CSVs · Drag the crosshair to set quadrant gates"
                :"Drop raw .fcs files · adjust the Cells & Singlets polygon gates · plot the gated singlet population"}
            </p>
          </div>
          <div style={{display:"flex",flexDirection:"column",alignItems:"flex-end",gap:8}}>
            {/* Mode switch */}
            <div style={{display:"flex",gap:4,padding:4,background:"#EEF0F3",borderRadius:10,border:"1px solid #E5E7EB"}}>
              {modeBtn("histogram",<HistIcon/>,"Histogram")}
              {modeBtn("quadrant",<QuadIcon/>,"Quadrant")}
              {modeBtn("analysis",<AnalysisIcon/>,"Analysis")}
            </div>
          </div>
        </div>
      </div>

      {/* Shared palette + gate-color + clear bar */}
      {samples.length>0&&mode!=="analysis"&&(
        <div style={{maxWidth:1200,margin:"0 auto 12px",background:"white",borderRadius:10,border:"1px solid #E5E7EB",padding:"10px 14px",display:"flex",alignItems:"center",gap:"10px 16px",flexWrap:"wrap"}}>
          <PalettePicker value={paletteName} onChange={applyPalette}/>
          <div style={{display:"flex",alignItems:"center",gap:8,marginLeft:"auto",flexWrap:"wrap"}}>
            <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em"}}>Gate</span>
            <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
              <button onClick={()=>setShowGate(true)} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:showGate?"#EFF6FF":"white",color:showGate?"#3B82F6":"#9CA3AF"}} title="Show the gate line on all plots">Show</button>
              <button onClick={()=>setShowGate(false)} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:!showGate?"#EFF6FF":"white",color:!showGate?"#3B82F6":"#9CA3AF"}} title="Hide the gate line on all plots">Hide</button>
            </div>
            <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em"}}>% labels</span>
            <div style={{display:"flex",borderRadius:6,border:"1px solid #E5E7EB",overflow:"hidden"}}>
              <button onClick={()=>setShowPct(true)} style={{padding:"4px 10px",border:"none",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:showPct?"#EFF6FF":"white",color:showPct?"#3B82F6":"#9CA3AF"}} title="Show percent labels on all plots">Show</button>
              <button onClick={()=>setShowPct(false)} style={{padding:"4px 10px",border:"none",borderLeft:"1px solid #E5E7EB",fontSize:11,fontWeight:600,cursor:"pointer",fontFamily:"var(--ff)",background:!showPct?"#EFF6FF":"white",color:!showPct?"#3B82F6":"#9CA3AF"}} title="Hide percent labels on all plots">Hide</button>
            </div>
            <span style={{fontSize:10,fontWeight:700,color:"#9CA3AF",textTransform:"uppercase",letterSpacing:"0.06em"}}>Color</span>
            <ColorPicker color={gateColor} onChange={setGateColor} align="right"/>
            <button onClick={()=>setMergeOpen(true)} title="Unify differently-named channels that are the same detector" style={{padding:"4px 12px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:12,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>Merge channels</button>
            <button onClick={clearAll} style={{padding:"4px 14px",borderRadius:6,border:"1px solid #E5E7EB",background:"white",fontSize:12,color:"#6B7280",cursor:"pointer",fontFamily:"var(--ff)"}}>Clear all</button>
          </div>
        </div>
      )}

      {mode==="analysis"
        ?<AnalysisMode gateColor={gateColor} showGate={showGate} showPct={showPct} onToggleGate={()=>setShowGate(g=>!g)} onTogglePct={()=>setShowPct(p=>!p)} paletteName={paletteName}/>
        :mode==="histogram"
        ?<HistogramMode samples={samples} allHeaders={allHeaders} colors={colors} updateSampleName={updateSampleName} updateColor={updateColor} removeSample={removeSample} swapSamples={swapSamples} moveSample={moveSample} hidden={hidden} toggleHidden={toggleHidden} uploadProps={uploadProps} gateColor={gateColor} showGate={showGate} onToggleGate={()=>setShowGate(g=>!g)} showPct={showPct} onTogglePct={()=>setShowPct(p=>!p)}/>
        :<QuadrantMode samples={samples} allHeaders={allHeaders} colors={colors} updateSampleName={updateSampleName} updateColor={updateColor} removeSample={removeSample} swapSamples={swapSamples} moveSample={moveSample} uploadProps={uploadProps} gateColor={gateColor} showGate={showGate} showPct={showPct}/>}

      {samples.length===0&&mode!=="analysis"&&<div style={{maxWidth:640,margin:"40px auto",color:"#6B7280",fontSize:13,lineHeight:1.75}}>
        <div style={{background:"white",border:"1px solid #E5E7EB",borderRadius:12,padding:"18px 22px"}}>
          <div style={{fontWeight:700,color:"#111827",marginBottom:8,fontSize:14}}>What to upload</div>
          <p style={{margin:"0 0 8px"}}>One CSV per sample, already <b>gated to the population you want to plot</b> (typically live singlets). Flume plots the events exactly as exported — it does not gate scatter for you in this view.</p>
          <p style={{margin:"0 0 8px"}}>In FlowJo: select the gated population → <b>File → Export / Concatenate → Export</b> → format <b>CSV — Channel Values</b>, and include {mode==="quadrant"?<>the <b>two fluorescence channels</b> you'll put on X and Y (plus a viability/scatter channel if you want to gate on it)</>:<>each <b>fluorescence channel</b> you want to plot</>}.</p>
          <p style={{margin:0,color:"#9CA3AF",fontSize:12}}>Have raw, ungated <b>.fcs</b> files instead? Use the <b>Analysis</b> tab to gate cells &amp; singlets first, then plot.</p>
        </div>
      </div>}

      {mergeOpen&&<MergeChannelsModal headers={allHeaders} onMerge={mergeChannels} onClose={()=>setMergeOpen(false)}/>}
    </div>
  );
}
