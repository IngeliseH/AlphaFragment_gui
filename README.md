# README for web_fragmenter

This is a standalone web app for interactive protein fragmentation using AlphaFragment.

## Setup

1. Create and activate a virtual environment:
   python3 -m venv env
   source env/bin/activate

2. Install requirements:
   pip install -r requirements.txt

3. Run the app:
   python app.py

## Usage
- Open the web page in your browser.
- Enter a UniProt ID to fragment and visualize interactively.

## Notes
- This app does not modify the main AlphaFragment codebase.
- Plots are interactive and rendered in-browser using Plotly.
