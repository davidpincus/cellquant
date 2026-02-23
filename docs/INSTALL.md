# Installation Guide

This guide assumes you have never used a terminal before. Every step is explicit. If something doesn't work, copy the error message and paste it to your AI assistant (Claude, ChatGPT, etc.) — it will know what to do.

## Step 1: Open your terminal

**Mac:** Press `Cmd + Space`, type "Terminal", press Enter. A window with a text prompt will appear. This is where you'll type commands.

**Windows:** Press `Win`, type "PowerShell", press Enter. (Note: some commands below may differ slightly on Windows. If something fails, tell your AI assistant you're on Windows and paste the error.)

**Linux:** You probably already know how to do this.

## Step 2: Install Miniforge (a Python distribution)

You need Python and a way to manage packages. Miniforge is the easiest option.

**Check if you already have conda:**
```bash
conda --version
```
If this prints a version number (e.g., `conda 24.x.x`), skip to Step 3.

**If you don't have conda, install Miniforge:**

Go to https://github.com/conda-forge/miniforge#download and download the installer for your operating system. Then:

**Mac/Linux:**
```bash
bash ~/Downloads/Miniforge3-$(uname)-$(uname -m).sh
```
Follow the prompts (press Enter to accept defaults, type "yes" when asked). Close and reopen your terminal after installation.

**Windows:**
Run the downloaded `.exe` installer. Accept all defaults.

**Verify it worked:**
```bash
conda --version
```
You should see a version number.

## Step 3: Create the cellquant environment

This creates an isolated Python environment with all the packages `cellquant` needs. This way, it won't interfere with anything else on your computer.

**If you have the `environment.yml` file** (from the GitHub repo):
```bash
conda env create -f environment.yml
```

**If you don't have it yet**, create the environment manually:
```bash
conda create -n cellquant python=3.11 -y
conda activate cellquant
pip install cellpose scikit-image numpy pandas matplotlib scipy pyyaml tifffile
```

## Step 4: Activate the environment

Every time you open a new terminal to run `cellquant`, you need to activate the environment first:

```bash
conda activate cellquant
```

Your terminal prompt should change to show `(cellquant)` at the beginning. If it shows `(base)` or nothing, the environment is not active.

## Step 5: Download cellquant.py

**Option A: Clone the whole repository** (if you have git):
```bash
git clone https://github.com/[username]/AVCGTIA.git
cd AVCGTIA
```

**Option B: Download just the script:**
Go to [GitHub URL], click on `cellquant.py`, click the "Raw" button, then `Cmd+S` (Mac) or `Ctrl+S` (Windows) to save the file. Remember where you saved it.

## Step 6: Verify the installation

```bash
conda activate cellquant
python cellquant.py --help
```

You should see a help message listing all available options. If you see an error about a missing package, install it:

```bash
pip install [package-name]
```

## Step 7: Download example data (optional)

To test the pipeline before using your own images:

```bash
cd AVCGTIA
ls example_data/mammalian_SGs/
```

You should see `.tif` files. If the example data folder is empty, download it from [Zenodo/FigShare link].

## Troubleshooting installation

**"conda: command not found"**
Close your terminal completely and reopen it. If it still doesn't work, Miniforge may not have been added to your PATH. Your AI assistant can help you fix this — tell it your operating system and paste the error.

**"ModuleNotFoundError: No module named 'cellpose'"**
You forgot to activate the environment. Run `conda activate cellquant` first.

**Cellpose is very slow to install**
This is normal — it installs PyTorch, which is a large download. Wait for it to finish (5–15 minutes on a typical connection).

**"ERROR: Could not find a version that satisfies the requirement..."**
Your Python version might be too old. Check with `python --version`. You need 3.11 or higher.

**Apple Silicon Mac: "BFloat16 is not supported on MPS"**
This is handled automatically by `cellquant.py`. It will detect your Mac and use CPU mode. No action needed. You'll see a warning message — this is expected and normal.

## You're done!

Proceed to [Tutorial 1: Mammalian Stress Granules](TUTORIAL_1_mammalian_SGs.md) to run your first analysis.

## Getting help

If you're stuck at any point:
1. Copy the exact error message from your terminal
2. Paste it to your AI assistant (Claude, ChatGPT, etc.)
3. Tell the AI: "I'm trying to install cellquant for fluorescence image analysis and got this error. I'm on [Mac/Windows/Linux]."

The AI will know what to do.
