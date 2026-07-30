import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// `base` must match the GitHub Pages repo subpath so built asset URLs resolve at
// https://<user>.github.io/flume/ . Override with BASE_PATH for a custom domain or root deploy.
export default defineConfig({
  base: process.env.BASE_PATH || "/flume/",
  plugins: [react()],
});
