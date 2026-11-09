# Giving this to someone else

There are two ways. Pick by whether the repository is public.

## Simplest: public repo, one file, double-click

The repository has to be **public** for this, because the installer downloads
from GitHub without an account.

Send them this link:

```
https://github.com/YoussefHarrabii/chem-tool/raw/main/chemdraft-setup.cmd
```

They save the file and double-click it. It installs everything and puts the
skill where opencode finds it. To update later, they double-click the same file
again; it always fetches the current version.

If the file opens in a browser instead of downloading, use right-click →
Save link as.

## Private repo: send the zip

If the repository stays private, the download needs an account, so do not use
the one-file installer. Instead:

1. Build the zip (`tools/make_zip.ps1` or zip the folder yourself).
2. Send it by email, WhatsApp or Drive.
3. They extract it anywhere and double-click `install.cmd` inside.

Updates are then manual: send a new zip.

They still need Python 3.10 or newer once. `install.cmd` will say so if it is
missing.

## For developers: git

If someone will pull updates, clone instead:

```
git clone https://github.com/YoussefHarrabii/chem-tool.git C:\Tools\chemdraft
cd C:\Tools\chemdraft
install.cmd
```

Private repositories need a one-time browser sign-in on first clone, or a
personal access token. After that `update.cmd` pulls and refreshes the skill.

## What each installer does

- installs RDKit, PyMuPDF and Pillow (RDKit is about 35 MB);
- installs this package in editable mode;
- copies `skills/chemdraft/SKILL.md` to
  `%USERPROFILE%\.config\opencode\skills\chemdraft\`;
- draws a test molecule and reports whether it worked.

## The prompt for opencode

Once installed, nothing else is needed: the skill is discovered automatically.
If you want to be explicit, paste this:

---

Use the `chemdraft` skill for every chemistry figure. Draw to a PNG under
`figures/` in the current project, tell me the formula the tool reported and the
file path, and pass on any warnings it returns. If the tool refuses, say so and
explain what it said; do not describe a figure that was not written.

---
