// Headless Chromium on the studio served by serve.jl (port 9384), driven through a command file: write one JSON command
// per file to bcmd.json, the result appears in bout.json. Commands:
//   {"op": "shot", "path": "x.png"}                 screenshot of the page
//   {"op": "click", "widget": "tab_files"}           real mouse click at the centre of a widget (bbox from state.json)
//   {"op": "clickxy", "x": 10, "y": 20}              click at figure coordinates (Makie: origin bottom left)
//   {"op": "clickpage", "x": 217, "y": 388}          click at page coordinates (as in a screenshot)
//   {"op": "type", "text": "abc", "enter": true}     keyboard input (after a click into a text box)
//   {"op": "key", "key": "Backspace", "n": 5}        press a key n times (Playwright key names)
//   {"op": "reload"}                                 reload the page (new WGLMakie session)
//   {"op": "quit"}
//     cd tools/studio && npm install && node browser.js
const { chromium } = require("playwright-core");
const fs = require("fs"), path = require("path");
const DIR = __dirname;
const sleep = (ms) => new Promise(r => setTimeout(r, ms));
const state = () => JSON.parse(fs.readFileSync(path.join(DIR, "state.json"), "utf8"));
const out = (o) => { fs.writeFileSync(path.join(DIR, "bout.json.tmp"), JSON.stringify(o)); fs.renameSync(path.join(DIR, "bout.json.tmp"), path.join(DIR, "bout.json")); };

(async () => {
  const browser = await chromium.launch({ executablePath: "/usr/bin/chromium",
    args: ["--no-sandbox", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"] });
  const page = await browser.newPage({ viewport: { width: 1900, height: 1000 } });
  const log = [];
  page.on("console", m => { if (m.type() === "error") log.push("[browser] " + m.text()); });
  page.on("pageerror", e => log.push("[pageerror] " + e.message));
  const open = async () => { await page.goto("http://localhost:9384/"); await page.waitForSelector("canvas", { timeout: 180000 }); await sleep(3000); };
  await open();
  async function toPage(x, y) {        // figure coordinates -> page coordinates
    const s = state();
    const r = await page.evaluate(() => { const b = document.querySelector("canvas").getBoundingClientRect(); return [b.left, b.top, b.height]; });
    const k = r[2] / s.fig[1];
    return [r[0] + x * k, r[1] + (s.fig[1] - y) * k];
  }
  out({ ready: true });
  const cmdFile = path.join(DIR, "bcmd.json");
  while (true) {
    if (fs.existsSync(cmdFile)) {
      const c = JSON.parse(fs.readFileSync(cmdFile, "utf8")); fs.unlinkSync(cmdFile);
      try {
        if (c.op === "shot") { await page.screenshot({ path: path.resolve(DIR, c.path), fullPage: !!c.full, timeout: c.timeout || 30000 }); out({ ok: c.path, log: log.splice(0) }); }
        else if (c.op === "click" || c.op === "clickxy") {
          let x = c.x, y = c.y;
          if (c.op === "click") { const b = state().widgets[c.widget]; if (!b) throw new Error("no widget " + c.widget); x = b[0] + b[2] / 2; y = b[1] + b[3] / 2; }
          const [px, py] = await toPage(x, y); await page.mouse.click(px, py); await sleep(c.wait || 500);
          out({ ok: c.op, at: [px, py], log: log.splice(0) });
        }
        else if (c.op === "clickpage") { await page.mouse.click(c.x, c.y); await sleep(c.wait || 500); out({ ok: "clickpage", log: log.splice(0) }); }
        else if (c.op === "key") { for (let i = 0; i < (c.n || 1); i++) await page.keyboard.press(c.key); await sleep(300); out({ ok: "key" }); }
        else if (c.op === "type") { await page.keyboard.type(c.text); if (c.enter) await page.keyboard.press("Enter"); await sleep(300); out({ ok: "type" }); }
        else if (c.op === "reload") { await open(); out({ ok: "reload" }); }
        else if (c.op === "quit") { out({ ok: "quit" }); await browser.close(); process.exit(0); }
        else throw new Error("unknown op " + c.op);
      } catch (e) { out({ error: e.message, log: log.splice(0) }); }
    }
    await sleep(100);
  }
})();
