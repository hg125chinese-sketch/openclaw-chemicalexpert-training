# Security & Privacy Guidelines

This repository is intended for public sharing.

## Do not commit secrets

Never commit any of the following:

- API keys, OAuth tokens, access tokens, refresh tokens
- Private keys / certificates
- Passwords, cookies, session files
- `.env` files containing credentials

If you need to show an API call, use **environment variables** (placeholders), e.g. `API_KEY=$YOUR_API_KEY`.

## Do not commit personal environment details

Avoid committing information that can identify a person or a private machine:

- Absolute paths to personal home directories (e.g., `/home/<user>/...`, `C:\\Users\\<name>\\...`)
- IP addresses or internal hostnames
- Cloud resource identifiers tied to personal accounts

Use placeholders instead:

- `<OPENCLAW_WORKSPACE>`
- `<PATH_TO_QMD>`
- `<REPO_URL>`

## Windows download artifacts

Do not commit Windows Mark-of-the-Web metadata files such as:

- `*:Zone.Identifier`

These are not part of the project and should be deleted.

## Reporting

If you discover a leaked secret in Git history, rotate the credential immediately and rewrite history if required.
