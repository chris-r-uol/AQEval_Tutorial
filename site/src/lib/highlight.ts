/**
 * A very small syntax highlighter for the tutorial's R and Python.
 *
 * The code shown is a dozen lines of function calls, so a full highlighting
 * library would be most of the page's JavaScript for four colours. This knows
 * comments, strings, numbers and the handful of words that matter.
 */

const KEYWORDS = new Set(['TRUE', 'FALSE', 'NULL', 'True', 'False', 'None', 'import', 'library', 'from', 'as']);

const TOKEN = /(#.*$)|("[^"]*"|'[^']*')|(\b\d+(?:\.\d+)?\b)|([A-Za-z_][A-Za-z0-9_.]*)(?=\s*\()|([A-Za-z_][A-Za-z0-9_]*)/gm;

const escape = (text: string) =>
	text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

/** Returns HTML with `<span class="tok-…">` around each recognised token. */
export function highlight(code: string): string {
	let html = '';
	let position = 0;
	for (const match of code.matchAll(TOKEN)) {
		html += escape(code.slice(position, match.index));
		const [text, comment, string, number, call, word] = match;
		if (comment) html += `<span class="tok-comment">${escape(text)}</span>`;
		else if (string) html += `<span class="tok-string">${escape(text)}</span>`;
		else if (number) html += `<span class="tok-number">${text}</span>`;
		else if (call && !KEYWORDS.has(text)) html += `<span class="tok-call">${escape(text)}</span>`;
		else if (KEYWORDS.has(word ?? call)) html += `<span class="tok-keyword">${text}</span>`;
		else html += escape(text);
		position = match.index + text.length;
	}
	return html + escape(code.slice(position));
}
