// ─── FCS parser (FCS 3.0 / 3.1, list mode; datatypes F, D, I) ───
// Validated against Attune NxT (FCS3.1) and MACSQuant (FCS3.0) exports.
export function parseFCS(arrayBuffer){
  const dv = new DataView(arrayBuffer);
  const asciiAt = (s,e)=>{let o="";for(let i=s;i<=e;i++)o+=String.fromCharCode(dv.getUint8(i));return o;};
  const version = asciiAt(0,5);
  if(!/^FCS/.test(version)) throw new Error("Not an FCS file");
  const num = (s,e)=>parseInt(asciiAt(s,e).trim(),10);

  let tStart=num(10,17), tEnd=num(18,25), dStart=num(26,33), dEnd=num(34,41);

  // TEXT segment: first char is the delimiter
  const delim = String.fromCharCode(dv.getUint8(tStart));
  const textRaw = asciiAt(tStart+1, tEnd);
  const parts = textRaw.split(delim);
  const kv = {};
  for(let i=0;i+1<parts.length;i+=2){ kv[parts[i].trim().toUpperCase()] = parts[i+1]; }

  if(!dStart){ dStart = parseInt(kv["$BEGINDATA"],10); dEnd = parseInt(kv["$ENDDATA"],10); }

  const par = parseInt(kv["$PAR"],10);
  const tot = parseInt(kv["$TOT"],10);
  const dtype = (kv["$DATATYPE"]||"F").trim().toUpperCase();
  const byteord = (kv["$BYTEORD"]||"1,2,3,4").trim();
  const little = byteord.startsWith("1");

  const names=[], longNames=[], bits=[], ranges=[];
  for(let p=1;p<=par;p++){
    names.push((kv["$P"+p+"N"]||("P"+p)).trim());
    longNames.push((kv["$P"+p+"S"]||"").trim());
    bits.push(parseInt(kv["$P"+p+"B"],10));
    ranges.push(parseFloat(kv["$P"+p+"R"]));
  }
  // Integer data is read as 8/16/32-bit words. Other bit widths (e.g. 10/12/14/24-bit) would
  // misalign the entire DATA segment, so reject them loudly instead of returning garbage.
  if(dtype==="I"){
    for(let p=0;p<par;p++){
      if(bits[p]!==8&&bits[p]!==16&&bits[p]!==32)
        throw new Error("Unsupported integer bit width $P"+(p+1)+"B="+bits[p]+" (Flume reads 8/16/32-bit integer, or 32/64-bit float, list-mode FCS). Re-export as FCS 3.1 float, or export channel values as CSV.");
    }
  }

  // De-duplicate channel names (some files repeat)
  const seen={};
  const headers = names.map(n=>{ if(seen[n]!==undefined){seen[n]++;return n+"."+seen[n];} seen[n]=0; return n; });

  const columns = {}; headers.forEach(h=>columns[h]=new Float32Array(tot));
  let off = dStart;
  for(let e=0;e<tot;e++){
    for(let p=0;p<par;p++){
      let v;
      if(dtype==="F"){ v=dv.getFloat32(off,little); off+=4; }
      else if(dtype==="D"){ v=dv.getFloat64(off,little); off+=8; }
      else { const b=bits[p]; if(b===16){v=dv.getUint16(off,little);off+=2;} else if(b===8){v=dv.getUint8(off);off+=1;} else {v=dv.getUint32(off,little);off+=4;} }
      columns[headers[p]][e]=v;
    }
  }
  const longMap={}, rangeMap={};
  headers.forEach((h,i)=>{ if(longNames[i]) longMap[h]=longNames[i]; rangeMap[h]=ranges[i]||0; });
  return { version, par, tot, headers, columns, longMap, ranges:rangeMap, cyt:(kv["$CYT"]||"").trim() };
}

// Detect scatter channels by common naming conventions
export function detectScatter(headers){
  const find = re => headers.find(h=>re.test(h));
  return {
    fscA: find(/^FSC-?A$/i) || find(/FSC.*A/i) || find(/FSC/i),
    sscA: find(/^SSC-?A$/i) || find(/SSC.*A/i) || find(/SSC/i),
    fscH: find(/^FSC-?H$/i) || find(/FSC.*H/i) || find(/FSC/i),
  };
}
