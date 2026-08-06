# Installation Guide

This guide assumes you have never used a terminal before. Every step is explicit. The whole process takes about 15–20 minutes and uses roughly 4 GB of disk space (mostly Python packages and the Cellpose segmentation model).

If something doesn't work at any step, copy the error message and paste it to your AI assistant (Claude, ChatGPT, etc.) — it will know what to do.

## Step 1: Open your terminal

**Mac:** Press `Cmd + Space`, type "Terminal", press Enter. A window with a text prompt will appear. This is where you'll type commands for the rest of this guide.

**Windows:** We recommend installing WSL (Windows Subsystem for Linux) first — it gives you a Linux terminal inside Windows and avoids compatibility issues. Open PowerShell as Administrator and run:
```
wsl --install
```
Restart your computer, then open "Ubuntu" from the Start menu. This is your terminal. All commands below will work as written.

If you'd rather not install WSL, the pipeline can work in PowerShell directly, but some steps may need adaptation. Ask your AI assistant for help if you hit issues.

**Linux:** You probably already know how to do this.

## Step 2: Install a Python package manager

You need a way to install Python and the libraries cellquant depends on. Pick one of the two paths below. Both work; **uv is recommended** because it's faster and `uv run cellquant.py` handles environment setup for you on every run.

### Path A (recommended): uv

uv is a single fast tool that manages Python versions and dependencies. Install it with one command:

**Mac / Linux / WSL:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows PowerShell (if not using WSL):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Close and reopen your terminal so the install registers. Verify:
```bash
uv --version
```

That's it — `uv` will fetch Python and all dependencies the first time you run cellquant (see Step 4, Path A). No environment activation needed in subsequent terminals.

### Path B: conda (Miniforge)

If you'd rather use conda (e.g., you already have it from another project), install Miniforge:

**Check if you already have conda:**
```bash
conda --version
```
If this prints a version number, skip to Step 3.

Otherwise, go to https://github.com/conda-forge/miniforge#download and download the installer for your system. Then run it:

**Mac (Apple Silicon — M1/M2/M3/M4):**
```bash
bash ~/Downloads/Miniforge3-MacOSX-arm64.sh
```

**Mac (Intel):**
```bash
bash ~/Downloads/Miniforge3-MacOSX-x86_64.sh
```

**Linux / WSL:**
```bash
bash ~/Downloads/Miniforge3-Linux-x86_64.sh
```

Follow the prompts: press Enter to scroll through the license, type "yes" to accept, press Enter to accept the default install location, and type "yes" when asked to initialize conda.

**Important:** Close your terminal completely and reopen it after installation. Conda won't work until you do this.

**Verify it worked:**
```bash
conda --version
```
You should see a version number. If you see "command not found," try closing and reopening your terminal one more time.

## Step 3: Download cellquant

You need the cellquant repository, which contains the pipeline script, example data, and all documentation.

**Option A: If you have git** (most Macs do):
```bash
git clone https://github.com/davidpincus/cellquant.git
cd cellquant
```

Not sure if you have git? Try running `git --version`. If it prints a version number, you're good. On Mac, it may prompt you to install developer tools — say yes and wait for that to finish, then try the clone again.

**Option B: Download as a zip:**
Go to https://github.com/davidpincus/cellquant, click the green "Code" button, click "Download ZIP." Unzip the downloaded file and open your terminal in that folder:
```bash
cd ~/Downloads/cellquant-main
```

Either way, you should now be inside the cellquant folder. Verify by running:
```bash
ls cellquant.py
```
If it prints `cellquant.py`, you're in the right place. If it says "No such file," you need to `cd` into the right directory.

> **Moving the folder:** You can move the `cellquant-main` folder anywhere on your computer (e.g., your home directory or a projects folder). Just make sure to use the correct path when running commands. For example, if you moved it to `~/cellquant-main/`, you'd run `python ~/cellquant-main/cellquant.py ...`

## Step 4: Install dependencies

Pick the path that matches what you installed in Step 2. Both work; only one is needed.

### Path A (recommended): uv

cellquant has its dependencies declared inline at the top of `cellquant.py` (a PEP 723 script header). `uv` reads that header and sets everything up automatically the first time you run the script:

```bash
uv run cellquant.py --help
```

On the first run, uv will:
1. Download a matching Python interpreter (~30 MB)
2. Resolve and install numpy/scipy/scikit-image/matplotlib/pandas/PyYAML/tifffile/Cellpose/opencv/torch into an isolated cache (~3 GB)
3. Print the cellquant help text

Subsequent invocations reuse the cache and start in a fraction of a second. There is no environment to "activate" — `uv run cellquant.py …` always uses the right dependencies. If you prefer to materialize a project-local venv:

```bash
uv venv
uv pip install -r requirements.txt
source .venv/bin/activate
```

### Path B: conda

```bash
conda env create -f environment.yml
```

This will take 5–15 minutes depending on your internet connection — it's downloading Python, PyTorch, Cellpose, and other scientific computing packages.

**If this fails**, try using `mamba` (a faster alternative that comes with Miniforge):
```bash
mamba env create -f environment.yml
```

**If that also fails**, you can create the environment manually as a last resort:
```bash
conda create -n cellquant python=3.11 -y
conda activate cellquant
pip install -r requirements.txt
```

### Path C: plain pip venv (fallback)

```bash
python3.11 -m venv cellquant_env
source cellquant_env/bin/activate       # Mac / Linux
pip install -r requirements.txt
```

> **Important package names:** If you install dependencies one by one, note that the YAML library is `pip install pyyaml` (not `yaml`), and scikit-image is `pip install scikit-image` (not `skimage`). These naming mismatches are a common gotcha.

## Step 5: Verify the installation

**uv path:**
```bash
uv run cellquant.py --help
```

**conda or pip-venv path:**
```bash
conda activate cellquant   # or: source cellquant_env/bin/activate
python cellquant.py --help
```

You should see a help message listing all available options. Notes:

- For conda/pip-venv paths, your terminal prompt should now show `(cellquant)` (or your venv name) at the beginning. If not, the environment isn't active — re-run the activate command.
- **For conda/pip-venv, every new terminal you open** requires you to re-activate the environment first.
- For `uv run`, there is no activate step — just call `uv run cellquant.py …` directly.

If you see an error about a missing package, install it (pip or `uv add`):
```bash
pip install [package-name]
```

## Note: Cellpose model download

The first time you process images, Cellpose will download its segmentation model (~500 MB). This is a one-time download that requires an internet connection. If you're on a slow connection, be patient — it will show a progress bar. The `--help` flag works offline, but actual image processing needs the model downloaded once.

## Step 6: GPU setup (optional, recommended for 3D)

cellquant auto-detects GPUs at startup and uses them when available. 2D analysis runs comfortably on CPU; 3D segmentation is roughly an order of magnitude faster on a discrete GPU and is the regime where setting one up matters most.

Whatever you do here, cellquant always falls back to CPU when no GPU is available, so this step is never blocking.

### Windows / WSL with NVIDIA CUDA

We recommend CUDA 12.1 to match the PyTorch wheels Cellpose pins against. From your activated environment (or replace `pip` with `uv pip` if using `uv venv`):

```bash
pip install --upgrade --index-url https://download.pytorch.org/whl/cu121 torch
```

Verify CUDA is visible to PyTorch:
```bash
python -c "import torch; print('CUDA:', torch.cuda.is_available()); print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"
```

If CUDA is available, cellquant will use it automatically. To force CPU, pass `--no-gpu`.

If `torch.cuda.is_available()` returns `False` despite an NVIDIA card, the most common cause is a driver mismatch: install the NVIDIA driver for CUDA 12.x from https://developer.nvidia.com/cuda-12-1-0-download-archive, reboot, and reinstall the CUDA-flavored torch wheel.

### Linux with NVIDIA CUDA

Same as Windows above — install the CUDA-flavored torch wheel:
```bash
pip install --upgrade --index-url https://download.pytorch.org/whl/cu121 torch
```

Then verify with the same `torch.cuda.is_available()` check. The Linux NVIDIA driver typically comes from your distro's package manager (`apt install nvidia-driver-535` or similar). Reboot after driver installation.

### Apple Silicon (M1/M2/M3/M4) — MPS

PyTorch's Metal (MPS) backend is available on Apple Silicon, but Cellpose's cpsam Transformer ops are not yet fully supported under MPS. cellquant detects this at startup and **falls back to CPU automatically** — you'll see:

```
[warn] MPS GPU not supported by cpsam Transformer; using CPU
```

This is expected and not a bug; it's the right call until upstream MPS support catches up. Performance on Apple Silicon CPU is acceptable for 2D MIP analysis but slow for 3D — for production 3D work, prefer a Linux or Windows machine with an NVIDIA GPU.

### No GPU (CPU only)

2D analysis runs in seconds-to-minutes per image on CPU. 3D segmentation on CPU is several minutes per stack; tolerable for a few stacks at a time, painful for a dataset of dozens. If you only do 2D analysis, you can skip this step entirely.

## Step 7: Test with example data (recommended)

Run a quick test to make sure everything works before using your own images:

```bash
python cellquant.py example_data/mammalian_SGs/ \
  "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
  --cell-type mammalian \
  --out example_data/mammalian_SGs/test_output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

This processes the example mammalian stress granule dataset. It should take 1–5 minutes on CPU. When it finishes, check the QC overlays:

```bash
open example_data/mammalian_SGs/test_output/qc/    # Mac
```

If you see images with cyan cell outlines and yellow nuclear outlines drawn over your fluorescence images, everything is working. Proceed to [Tutorial 1](TUTORIAL_1_mammalian_SGs.md).

## Troubleshooting installation

**"conda: command not found"**
Close your terminal completely and reopen it. If it still doesn't work, Miniforge wasn't added to your PATH. Tell your AI assistant your operating system and paste the error — it can walk you through fixing this.

**"ModuleNotFoundError: No module named 'cellpose'"**
You forgot to activate the environment. Run `conda activate cellquant` first.

**Cellpose is very slow to install**
Normal — it's downloading PyTorch, which is about 2 GB. Wait for it to finish.

**"ERROR: Could not find a version that satisfies the requirement..."**
Your Python version might be too old. Check with `python --version`. You need 3.11 or higher.

**Apple Silicon Mac: "BFloat16 is not supported on MPS" or "MPS GPU not supported"**
This is handled automatically. The pipeline detects your Mac and uses CPU mode. You'll see a warning message — this is expected and everything will work correctly, just a bit slower.

**"RuntimeError: Numpy is not available" / "numpy.core.multiarray failed to import" / numpy 2.x errors**
Numpy 2.0+ is incompatible with the versions of PyTorch and Cellpose that cellquant uses. This happens when `pip` upgrades numpy to 2.x during installation. Fix it by pinning numpy and opencv:
```bash
pip install "numpy>=1.24,<2.0" "opencv-python-headless<4.10"
```
> **zsh users (default Mac shell):** The quotes around the version specs are required — without them, `<` is interpreted as a shell redirect and the command will fail silently or create junk files.

**"git: command not found" (Mac)**
Run `xcode-select --install` to install developer tools, then try again.

## You're done!

Proceed to [Tutorial 1: Mammalian Stress Granules](TUTORIAL_1_mammalian_SGs.md) to learn how the pipeline works, or jump to the [Quick Start](QUICKSTART.md) if you just want to run it.

## Getting help

If you're stuck at any point:
1. Copy the exact error message from your terminal
2. Paste it to your AI assistant (Claude, ChatGPT, etc.)
3. Tell the AI: "I'm trying to install cellquant for fluorescence image analysis and got this error. I'm on [Mac/Windows/Linux]."

The AI will know what to do.
