# Giving this to someone else

They already have opencode, so the agent installs it. One message to paste, no
downloads, no double-clicking, no Windows warnings.

**First:** the repository must be public, otherwise the clone asks for a GitHub
account. One setting, on GitHub: Settings → General → Danger Zone → Change
visibility → Public. Nothing in this repository is sensitive.

**Then send them this, and nothing else:**

---

Install chemdraft so you can draw chemistry diagrams for me:

1. `git clone https://github.com/YoussefHarrabii/chem-tool.git %USERPROFILE%\Tools\chemdraft`
   If the clone asks to sign in, tell me and wait for me to finish that.
2. In that folder run `install.cmd`.
3. Check it worked:
   `chemdraft smiles "CCO" -o "%TEMP%\check.png"`
   It must print JSON containing `"ok": true`.
```
4. Tell me the version it printed, then tell me to start a new chat so the
   skill is picked up.
```

From then on, use the chemdraft skill for every chemistry figure I ask for.
Update it when I say so with `update.cmd` in the same folder.

---

That is the whole setup. The agent clones, installs, registers the skill, and
verifies. If Python is missing, it will say so and can install it.

Three things to expect on their side:

- They paste that into the **opencode chat**, not a terminal window.
- opencode may ask permission to run commands. Allow it, once.
- **A skill loads when a session starts**, so the install only takes effect in a
  new chat. That is why step 4 says to start one.

## Updating

They say "update chemdraft". The agent runs `update.cmd`, which pulls the latest
version and refreshes the skill. No action from them beyond the message.

## Why not a double-click installer

A `.cmd` file downloaded from the internet gets flagged by SmartScreen and
Windows Defender, which is worse for a non-technical person than one paste. Files
obtained by `git clone` carry no such mark, so `install.cmd` runs cleanly when
the agent uses it.

## If the repository must stay private

Then the clone needs a sign-in on first use. Add them under Settings →
Collaborators, and expect one browser sign-in during step 1. Everything after
that is identical.

## What gets installed

- RDKit, PyMuPDF, Pillow (RDKit is about 35 MB);
- this package in editable mode, so `git pull` updates it;
- `skills/chemdraft/SKILL.md` copied to
  `%USERPROFILE%\.config\opencode\skills\chemdraft\`;
- a test figure written to `%TEMP%`.
