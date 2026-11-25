"""Every message the user typed in a Claude Code session transcript, in order, with timestamps and compaction marks.

Part of the catch-up the user asked for after every context compaction (project rule, 2026-11-25; see CLAUDE.md):
re-read all user messages, including corrections and decisions, before continuing. Tool results, system reminders,
task notifications and slash-command echoes are left out; each message is cut at --max-chars.

Usage: python extract_user_messages.py <transcript.jsonl> [--since-last-compaction] [--max-chars 1500]
       (transcripts: ~/.claude/projects/<project>/<session-id>.jsonl)
Output: printed, one line per message; "--- COMPACTION at <time> ---" where the context was compacted
"""
import argparse
import json


def messages(path, max_chars):
    out = []
    for line in open(path):
        try:
            d = json.loads(line)
        except ValueError:
            continue
        if d.get("isCompactSummary"):
            out.append(f"--- COMPACTION at {d.get('timestamp', '')[:16]} ---")
            continue
        if d.get("type") != "user" or d.get("isMeta"):
            continue
        c = d.get("message", {}).get("content")
        if isinstance(c, list):
            txt = " ".join(x.get("text", "") for x in c if isinstance(x, dict) and x.get("type") == "text")
        else:
            txt = c or ""
        txt = txt.strip()
        if (not txt or txt.startswith("<system-reminder>") or "<task-notification>" in txt
                or txt.startswith("<local-command") or txt.startswith("<command-name>")):
            continue
        out.append(f"[{d.get('timestamp', '')[:16]}] {txt[:max_chars]}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("transcript")
    ap.add_argument("--since-last-compaction", action="store_true")
    ap.add_argument("--max-chars", type=int, default=1500)
    args = ap.parse_args()
    out = messages(args.transcript, args.max_chars)
    if args.since_last_compaction:
        marks = [i for i, s in enumerate(out) if s.startswith("--- COMPACTION")]
        if marks:
            out = out[marks[-1]:]
    print("\n".join(out))
    print(f"\n({sum(not s.startswith('---') for s in out)} user messages, "
          f"{sum(s.startswith('---') for s in out)} compactions shown)")


if __name__ == "__main__":
    main()
