// =============================================================================
// Disjunctive Cuts Lab — desktop shell.
//
// Starts the compiled Julia server (build/server, produced by
// build/build_server.jl) on a free local port, shows a splash screen until it
// answers, then loads the web interface in a native window. The server is
// stopped when the app quits.
//
// Development without a compiled server: set CUTS_DEV_JULIA=1 to run
// `julia --project=app app/server.jl` instead (Julia must be on PATH).
// =============================================================================

const { app, BrowserWindow, dialog, shell } = require("electron");
const { spawn } = require("child_process");
const fs = require("fs");
const http = require("http");
const net = require("net");
const path = require("path");

const REPO_ROOT = path.join(__dirname, "..");
const EXE = process.platform === "win32" ? "cuts-server.exe" : "cuts-server";
const STARTUP_TIMEOUT_MS = 180_000;

let server = null;          // child process running the Julia server
let serverLog = "";         // last output of the server, shown if it fails
let quitting = false;

// ---- Locating resources -------------------------------------------------------

/** Path of a bundled resource: inside the installed app, or in the repo while developing. */
function resource(...parts) {
  return app.isPackaged ? path.join(process.resourcesPath, ...parts) : null;
}

function serverCommand(port) {
  // CUTS_WATCH_STDIN: the server exits when our end of its stdin pipe closes,
  // so it stops even if this app crashes or is force-quit.
  const env = { ...process.env, PORT: String(port), CUTS_PLOTLY_JS: plotlyPath(), CUTS_WATCH_STDIN: "1" };
  if (process.env.CUTS_DEV_JULIA) {
    return { cmd: "julia", args: ["--project=app", "app/server.jl"], cwd: REPO_ROOT, env };
  }
  const bin = resource("server", "bin", EXE) ?? path.join(REPO_ROOT, "build", "server", "bin", EXE);
  if (!fs.existsSync(bin)) {
    throw new Error(`Compiled server not found at\n${bin}\n\nBuild it with:\n  julia --project=build build/build_server.jl\n` +
                    `or run with CUTS_DEV_JULIA=1 to use a local Julia installation.`);
  }
  return { cmd: bin, args: [], cwd: path.dirname(bin), env };
}

function plotlyPath() {
  return resource("plotly.min.js") ?? require.resolve("plotly.js-dist-min/plotly.min.js");
}

// ---- Server lifecycle ---------------------------------------------------------

/** Ask the OS for a free TCP port on localhost. */
function freePort() {
  return new Promise((resolve, reject) => {
    const srv = net.createServer();
    srv.unref();
    srv.on("error", reject);
    srv.listen(0, "127.0.0.1", () => {
      const { port } = srv.address();
      srv.close(() => resolve(port));
    });
  });
}

function startServer(port) {
  const { cmd, args, cwd, env } = serverCommand(port);
  server = spawn(cmd, args, { cwd, env, windowsHide: true });
  const keep = (chunk) => { serverLog = (serverLog + chunk.toString()).slice(-8000); };
  server.stdout.on("data", keep);
  server.stderr.on("data", keep);
  return server;
}

/** Resolve once GET /health answers 200; reject if the server exits or times out. */
function waitForServer(port, child) {
  const deadline = Date.now() + STARTUP_TIMEOUT_MS;
  return new Promise((resolve, reject) => {
    let exited = false;
    child.once("exit", (code) => { exited = true; reject(new Error(`The solver stopped during startup (exit code ${code}).`)); });
    const poll = () => {
      if (exited) return;
      if (Date.now() > deadline) return reject(new Error("The solver did not start within 3 minutes."));
      const req = http.get({ host: "127.0.0.1", port, path: "/health", timeout: 1000 }, (res) => {
        res.resume();
        res.statusCode === 200 ? resolve() : setTimeout(poll, 300);
      });
      req.on("error", () => setTimeout(poll, 300));
      req.on("timeout", () => { req.destroy(); setTimeout(poll, 300); });
    };
    poll();
  });
}

function stopServer() {
  if (server && server.exitCode === null) server.kill();
  server = null;
}

// ---- Window -------------------------------------------------------------------

async function createWindow() {
  const win = new BrowserWindow({
    width: 1360, height: 900, minWidth: 720, minHeight: 560,
    title: "Disjunctive Cuts Lab", backgroundColor: "#f6f7f9", show: false,
    webPreferences: { contextIsolation: true, sandbox: true },
  });
  win.removeMenu?.();
  win.once("ready-to-show", () => win.show());

  // Links to other sites (e.g. documentation) open in the user's browser.
  win.webContents.setWindowOpenHandler(({ url }) => { shell.openExternal(url); return { action: "deny" }; });

  await win.loadFile(path.join(__dirname, "loading.html"));

  try {
    const port = await freePort();
    const child = startServer(port);
    await waitForServer(port, child);
    child.on("exit", (code) => {
      if (quitting) return;
      dialog.showErrorBox("Solver stopped", `The solver process exited unexpectedly (code ${code}).\n\n${serverLog.slice(-2000)}`);
      app.quit();
    });
    await win.loadURL(`http://127.0.0.1:${port}/`);
  } catch (err) {
    dialog.showErrorBox("Could not start the solver", `${err.message}\n\n${serverLog.slice(-2000)}`);
    app.quit();
  }
}

// One instance at a time: a second launch focuses the existing window.
if (!app.requestSingleInstanceLock()) {
  app.quit();
} else {
  app.on("second-instance", () => {
    const [win] = BrowserWindow.getAllWindows();
    if (win) { if (win.isMinimized()) win.restore(); win.focus(); }
  });
  app.whenReady().then(createWindow);
  app.on("window-all-closed", () => app.quit());
  app.on("before-quit", () => { quitting = true; stopServer(); });
}
