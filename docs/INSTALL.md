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

## Step 2: Install Miniforge (a Python distribution)

You need Python and a package manager called `conda`. Miniforge is the easiest way to get both.

**Check if you already have conda:**
```bash
conda --version
```
If this prints a version number (e.g., `conda 24.x.x`), skip to Step 3.

**If you don't have conda, install Miniforge:**

Go to https://github.com/conda-forge/miniforge#download and download the installer for your system. Then run it:

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
git clone https://github.com/pincus-lab/cellquant.git
cd cellquant
```

Not sure if you have git? Try running `git --version`. If it prints a version number, you're good. On Mac, it may prompt you to install developer tools — say yes and wait for that to finish, then try the clone again.

**Option B: Download as a zip:**
Go to https://github.com/pincus-lab/cellquant, click the green "Code" button, click "Download ZIP." Unzip the downloaded file and open your terminal in that folder:
```bash
cd ~/Downloads/cellquant-main
```

Either way, you should now be inside the cellquant folder. Verify by running:
```bash
ls cellquant.py
```
If it prints `cellquant.py`, you're in the right place. If it says "No such file," you need to `cd` into the right directory.

## Step 4: Create the cellquant environment

This creates an isolated Python environment with all the packages cellquant needs. It won't interfere with anything else on your computer.

```bash
conda env create -f environment.yml
```

This will take 5–15 minutes depending on your internet connection — it's downloading Python, PyTorch, Cellpose, and other scientific computing packages. Let it finish.

**If this fails**, try using `mamba` (a faster alternative that comes with Miniforge):
```bash
mamba env create -f environment.yml
```

**If that also fails**, you can create the environment manually as a last resort:
```bash
conda create -n cellquant python=3.11 -y
conda activate cellquant
pip install cellpose scikit-image numpy pandas matplotlib scipy pyyaml tifffile
```

## Step 5: Verify the installation

Activate the environment and check that cellquant runs:

```bash
conda activate cellquant
python cellquant.py --help
```

You should see a help message listing all available options. Two things to note:

- Your terminal prompt should now show `(cellquant)` at the beginning. If it shows `(base)` or nothing, the environment isn't active — run `conda activate cellquant` again.
- **Every time you open a new terminal** to use cellquant, you need to run `conda activate cellquant` first. This is easy to forget.

If you see an error about a missing package, install it:
```bash
pip install [package-name]
```

## Step 6: Test with example data (recommended)

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
