/**
 * Write the R script and the Python notebook from the tutorial's steps.
 *
 * The steps live in `src/lib/code.js`, which also feeds the website's code
 * panel, so the three always show the same code. Run after changing a step:
 *
 *     npm run examples            # rewrite the files
 *     npm run examples -- --check # fail if they are out of date (used in CI)
 */
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { DEFAULTS, plain, steps } from '../src/lib/code.js';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..', '..');
const SITE_URL = 'https://chris-r-uol.github.io/AQEval_Tutorial/';

/** Suggestions for what to change next. R's search for break points slows
 *  sharply as h falls, so the R version says how to keep it quick. */
const thingsToTry = (language) => [
	'Change "no2" to "nox". Do the changes fall on the same dates?',
	'Change the site from "bdma" to "led6" (Leeds Headingley Kerbside), a roadside site with no Clean Air Zone.',
	'Make the analysis more sensitive: change h from 0.3 to 0.2, then to 0.12, the value normally used for 8-hour data. Smaller values find more changes.' +
		(language === 'r'
			? ' In R the search then takes tens of minutes on 8-hour data, so change "8 hour" to "day" first: on daily data it takes a minute or two.'
			: ' On 8-hour data, measuring the changes then takes several minutes, so change "8 hour" to "day" first: on daily data it takes under a minute.'),
	'Use a background site as a control: the website shows the extra lines when you choose one under "Also account for".'
];

/** Break a paragraph into lines no longer than `width`. */
function wrap(text, width) {
	const lines = [];
	let line = '';
	for (const word of text.split(/\s+/)) {
		if (line && line.length + 1 + word.length > width) {
			lines.push(line);
			line = word;
		} else {
			line = line ? `${line} ${word}` : word;
		}
	}
	if (line) lines.push(line);
	return lines;
}

const comment = (text) => wrap(text, 76).map((line) => `# ${line}`);

function rScript() {
	const out = [
		"# Finding changes in Bradford's air quality with AQEval",
		'#',
		...comment(`This is the R version of the interactive tutorial at ${SITE_URL}`),
		'#',
		'# How to use it',
		...comment(
			'Click on the first line of code, then press Ctrl+Enter (Cmd+Enter on a Mac) to run it and move to the next. ' +
				'Plots appear in the Plots pane and tables in the Console.'
		),
		'#',
		...comment(
			'This file is generated from site/src/lib/code.js, so that it always matches the website. ' +
				'Change it as much as you like while you work.'
		),
		''
	];
	steps(DEFAULTS, 'r').forEach((step, index) => {
		const heading = `# ${index + 1}. ${step.title} `;
		out.push(heading + '-'.repeat(Math.max(4, 76 - heading.length)));
		out.push(...comment(step.text));
		if (step.wait) out.push(...comment(step.wait));
		out.push(plain(step.code), '');
	});
	out.push('# Things to try ' + '-'.repeat(60));
	thingsToTry('r').forEach((item, index) => {
		const [first, ...rest] = wrap(item, 72);
		out.push(`# ${index + 1}. ${first}`, ...rest.map((line) => `#    ${line}`));
	});
	return out.join('\n') + '\n';
}

function notebook() {
	let count = 0;
	const cell = (type, source) => ({
		cell_type: type,
		id: `cell-${String(++count).padStart(2, '0')}`,
		metadata: {},
		// Jupyter stores a cell as a list of lines, each but the last ending in a newline.
		source: source.split('\n').map((line, index, all) => (index < all.length - 1 ? `${line}\n` : line)),
		...(type === 'code' ? { execution_count: null, outputs: [] } : {})
	});

	const cells = [
		cell(
			'markdown',
			[
				"# Finding changes in Bradford's air quality with AQEval",
				'',
				`This is the Python version of the [interactive tutorial](${SITE_URL}).`,
				'',
				'**How to use it.** Click on a grey code box and press **Shift+Enter** to run it and move to the next. ' +
					'Run the boxes in order, from the top. The first time, choose the **Python 3** kernel if you are asked.'
			].join('\n')
		)
	];
	steps(DEFAULTS, 'python').forEach((step, index) => {
		cells.push(
			cell('markdown', `## ${index + 1}. ${step.title}\n\n${step.text}${step.wait ? `\n\n*${step.wait}*` : ''}`)
		);
		cells.push(cell('code', plain(step.code)));
	});
	cells.push(
		cell(
			'markdown',
			['## Things to try', '', ...thingsToTry('python').map((item, index) => `${index + 1}. ${item}`)].join('\n')
		)
	);

	return (
		JSON.stringify(
			{
				cells,
				metadata: {
					kernelspec: { display_name: 'Python 3', language: 'python', name: 'python3' },
					language_info: { name: 'python' }
				},
				nbformat: 4,
				nbformat_minor: 5
			},
			null,
			1
		) + '\n'
	);
}

const files = {
	'R/aqeval_tutorial.R': rScript(),
	'python/aqeval_tutorial.ipynb': notebook()
};

const check = process.argv.includes('--check');
let stale = false;
for (const [path, content] of Object.entries(files)) {
	const target = resolve(ROOT, path);
	if (check) {
		if (!existsSync(target) || readFileSync(target, 'utf8') !== content) {
			console.error(`${path} is out of date: run "npm run examples" in site/`);
			stale = true;
		}
	} else {
		mkdirSync(dirname(target), { recursive: true });
		writeFileSync(target, content);
		console.log(`wrote ${path}`);
	}
}
process.exit(stale ? 1 : 0);
