# Setup prompt

Paste this into opencode on a fresh machine. The agent does the install; you do
the sign-in if the repository is private.

---

Set up and use the **chemdraft** chemistry drawing tool.

Repository: `https://github.com/YoussefHarrabii/chem-tool`

1. Check `git` and `python` are on PATH, and that Python is 3.10 or newer. If
   either is missing, tell me exactly what to install and stop.
2. Clone the repository to `%USERPROFILE%\Tools\chemdraft` if it is not there
   already. If the clone asks for credentials, tell me to finish the GitHub
   sign-in in the browser window, then continue.
3. Run `install.cmd` from that folder. On macOS or Linux, run
   `python -m pip install -e .` and copy `skills/chemdraft/SKILL.md` to
   `~/.config/opencode/skills/chemdraft/SKILL.md`.
4. Verify with:
   `chemdraft smiles "CCO" -o "%TEMP%\chemdraft-check.png"`
   It must print JSON containing `"ok": true`. If it does not, show me the error
   and stop.
5. Read `skills/chemdraft/SKILL.md` and follow it for every chemistry figure from
   now on.

From then on, whenever I ask for a molecule, a structure, a reaction or a
mechanism:

- produce a PNG with `chemdraft`, saved under `figures/` in the current project
  with a descriptive kebab-case name;
- tell me the formula it reported and the file path;
- pass on any warnings or `caveats` the tool returned;
- if the tool refuses, say so and explain what it said. Do not describe a figure
  that was not written, and do not retype a SMILES until a check passes.

To update later: run `update.cmd` in the same folder.

---

## If you want it private

The repository can be private. Then:

- add each person under **Settings → Collaborators**;
- the first clone opens a browser sign-in, once per machine;
- if the sign-in does not appear, either install **GitHub Desktop** and clone
  through it, or create a personal access token (Settings → Developer settings →
  Tokens) and use it as the password when git asks.

An agent cannot complete that sign-in for you, because it needs your account.

## Manual fallback

```
git clone https://github.com/YoussefHarrabii/chem-tool.git C:\Tools\chemdraft
cd C:\Tools\chemdraft
install.cmd
```

`install.cmd` installs RDKit, PyMuPDF and Pillow, installs the package in
editable mode, and copies the skill to
`%USERPROFILE%\.config\opencode\skills\chemdraft\`. Because the install is
editable, `git pull` updates the code; `update.cmd` also refreshes the skill
copy.
