# Git hooks

These hooks are versioned with the repo (unlike `.git/hooks/`, which is not).
Enable them in each clone with:

```bash
git config core.hooksPath .githooks
```

## pre-commit

Blocks commits that contain:

- Files under `output/` or `data/input_data/` (personal genetic data)
- `.env*` files other than `.env.example` (real credentials)
- `.dev.vars` (Wrangler secrets)
- Common credential patterns in staged content (`R2_SECRET_ACCESS_KEY=...`, `aws_secret_access_key=...`, `api_key=...`, etc.)
- Files larger than 10 MB

To bypass in a genuine emergency:

```bash
git commit --no-verify
```

Do not use `--no-verify` to silence a true positive — fix the underlying issue.
