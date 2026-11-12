# Contributing to KENRYU

Thank you for your interest in contributing to KENRYU. This project is maintained by a single author, but external contributions, bug reports, and feature requests are welcome and are part of how the project grows.

## How to report a bug

Open an issue at https://github.com/abrangel/Kenryu/issues and include:

1. **Steps to reproduce** — the miRNA panel, consensus mode, and year window used.
2. **Expected vs. actual behavior.**
3. **Logs** — browser console output for frontend errors, Hugging Face Space logs for backend errors.
4. **Environment** — browser and version, or Python/Docker version for local runs.

## How to suggest a feature

Open an issue with the `enhancement` label describing the use case: what research question it would help answer and what the expected output looks like. Roadmap items (isomiRs, DisGeNET, structured JSON export, panel comparison, heatmap, batch CSV mode, TargetScan v9.0) already have open issues — comment there instead of opening a duplicate.

## How to contribute code

1. **Fork** the repository.
2. Create a branch: `git checkout -b feat/feature-name` (or `fix/bug-name`).
3. Commit with descriptive messages following the project conventions: `feat:` · `fix:` · `docs:` · `refactor:` · `chore:`.
4. Push: `git push origin feat/feature-name`.
5. Open a **Pull Request** against `main`, describing what changed and how you tested it.

## Development setup

```bash
git clone https://github.com/abrangel/Kenryu.git
cd Kenryu
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
uvicorn kenryu_engine:app --host 0.0.0.0 --port 7860
```

## Support

- **Bug reports and feature requests:** GitHub Issues.
- **Questions about the science:** open a Discussion or an issue with the `question` label.
- **Security issues:** see `SECURITY.md` — do not report vulnerabilities in public issues.

## License

By contributing, you agree that your contributions will be distributed under the project's Free Academic License for Research and Education.
