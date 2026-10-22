# Put the demo online as a Hugging Face Space (about 10 minutes, optional)

1. Push this repository to GitHub first (the Dockerfile clones it from github.com/tammy357/ts-second-opinion).
2. On huggingface.co: New → Space → name `ts-second-opinion`, SDK **Docker**, licence MIT, hardware **CPU basic** (free).
3. Upload `README.md` and `Dockerfile` from this folder to the Space (Files → Add file → Upload files).
4. Wait for the build (a few minutes), then open the Space's URL.
5. Add a badge to the top of the main README (replace YOUR-NAME):
   `[![Try it on Hugging Face](https://huggingface.co/datasets/huggingface/badges/resolve/main/open-in-hf-spaces-sm.svg)](https://huggingface.co/spaces/YOUR-NAME/ts-second-opinion)`
   and put the link in your submission form and video description.

The exact command the Space runs (`marimo run apps/demo.py --host 0.0.0.0 --port 7860 --headless` after a fresh
`pip install -e .`) was tested from a clean clone. The Docker build itself has not been run on Hugging Face yet: if it
fails, the log on the Space page says why. The free CPU Space sleeps when idle and takes about a minute to wake.
