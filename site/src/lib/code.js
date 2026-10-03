/**
 * The tutorial's code, written once.
 *
 * The website's code panel, the R script students run in RStudio and the
 * Python notebook are all generated from the steps below, so the three cannot
 * drift apart. `scripts/build-examples.mjs` imports this file from Node to
 * write the script and the notebook, which is why it is plain JavaScript with
 * no browser or Svelte imports.
 *
 * Values that come from a dial on the website are wrapped with `dial()`. The
 * code panel turns those into highlights; everything else strips them.
 */

/** @typedef {'r' | 'python'} Language */

/**
 * @typedef {object} Settings
 * @property {string} site        AURN site code, e.g. "BDMA"
 * @property {string} pollutant   "no2" or "nox"
 * @property {boolean} deseason
 * @property {boolean} deweather
 * @property {string | null} background  null, "air_temp", or a control site's code
 * @property {string} averaging   openair-style period, e.g. "8 hour"
 * @property {number} h
 */

/**
 * @typedef {object} Step
 * @property {string} id
 * @property {string} title
 * @property {string} text   What the step does, in plain English.
 * @property {string} code   Source, with dial markers.
 */

export const YEARS = [2018, 2025];

/**
 * Where the dials, the R script and the notebook all start.
 *
 * h starts high on purpose. 8-hour data would normally be analysed with
 * h = 0.12, but starting at 0.3 gives a simple first result, lets students
 * find out for themselves what lowering h does, and keeps the last step to
 * about a minute in R rather than several.
 * @type {Settings}
 */
export const DEFAULTS = {
	site: 'BDMA',
	pollutant: 'no2',
	deseason: true,
	deweather: true,
	background: 'air_temp',
	averaging: '8 hour',
	h: 0.3
};

// Private-use characters, so a marker can never collide with real code.
const OPEN = '';
const SEPARATOR = '';
const CLOSE = '';

/** Wrap a value that is set by a dial. */
function dial(/** @type {keyof Settings} */ name, /** @type {string} */ text) {
	return `${OPEN}${name}${SEPARATOR}${text}${CLOSE}`;
}

/** Remove the dial markers, leaving runnable source. */
export function plain(/** @type {string} */ code) {
	return code.replace(new RegExp(`${OPEN}[^${SEPARATOR}]*${SEPARATOR}|${CLOSE}`, 'g'), '');
}

/**
 * Split marked-up code into runs of ordinary text and dial values.
 * @returns {Array<{ text: string, dial: string | null }>}
 */
export function segments(/** @type {string} */ code) {
	const pattern = new RegExp(`${OPEN}([^${SEPARATOR}]*)${SEPARATOR}([^${CLOSE}]*)${CLOSE}`, 'g');
	const parts = [];
	let position = 0;
	for (const match of code.matchAll(pattern)) {
		if (match.index > position) parts.push({ text: code.slice(position, match.index), dial: null });
		parts.push({ text: match[2], dial: match[1] });
		position = match.index + match[0].length;
	}
	if (position < code.length) parts.push({ text: code.slice(position), dial: null });
	return parts;
}

/** Whether the settings ask for any signal isolation at all. */
export function isolates(/** @type {Settings} */ settings) {
	return settings.deseason || settings.deweather || settings.background !== null;
}

/** Whether the background is another monitoring site rather than a column. */
export function usesControlSite(/** @type {Settings} */ settings) {
	return settings.background !== null && settings.background !== 'air_temp';
}

const bool = {
	r: (/** @type {boolean} */ value) => (value ? 'TRUE' : 'FALSE'),
	python: (/** @type {boolean} */ value) => (value ? 'True' : 'False')
};

/**
 * The tutorial as a list of steps, for one language and one set of dials.
 * @param {Settings} settings
 * @param {Language} language
 * @returns {Step[]}
 */
export function steps(settings, language) {
	const r = language === 'r';
	const site = dial('site', `"${settings.site.toLowerCase()}"`);
	const pollutant = dial('pollutant', `"${settings.pollutant}"`);
	const control = usesControlSite(settings);
	const isolating = isolates(settings);
	// The column the break analysis runs on.
	const target = isolating ? '"isolated"' : pollutant;
	const [first, last] = YEARS;

	/** @type {Step[]} */
	const list = [];

	list.push({
		id: 'packages',
		title: 'Load the packages',
		text: r
			? 'openair downloads UK air quality data and AQEval finds the changes in it.' +
				(control ? ' dplyr joins two tables together.' : '')
			: 'aqeval is the Python version of the R package AQEval. It downloads UK air quality data and finds the changes in it.',
		code: r
			? ['library(openair)', 'library(AQEval)', ...(control ? [dial('background', 'library(dplyr)')] : [])].join('\n')
			: 'import aqeval'
	});

	list.push({
		id: 'data',
		title: 'Get the data',
		text:
			`Download every hourly measurement from ${first} to ${last} for one monitoring site on the national network (AURN). ` +
			'Each site has a short code. The download also brings the wind speed, wind direction and air temperature for each hour.',
		code: r
			? `data <- importAURN(site = ${site}, year = ${first}:${last})`
			: `data = aqeval.download_aurn_data(\n    ${site}, ${first}, ${last}, source="aurn"\n)`
	});

	if (control) {
		const controlSite = dial('background', `"${/** @type {string} */ (settings.background).toLowerCase()}"`);
		list.push({
			id: 'control',
			title: 'Add a background site',
			text:
				'Download the same pollutant from a background site, away from busy roads, and line it up hour by hour with the first site. ' +
				'It shows what the air was doing across the whole region, whatever happened on this one road.',
			code: r
				? [
						`background <- importAURN(site = ${controlSite}, year = ${first}:${last})`,
						`background <- select(background, date, background = ${dial('pollutant', settings.pollutant)})`,
						'data <- left_join(data, background, by = "date")'
					].join('\n')
				: [
						`background = aqeval.download_aurn_data(\n    ${controlSite}, ${first}, ${last}, source="aurn"\n)`,
						`background = background[["date", ${pollutant}]].rename(\n    columns={${pollutant}: "background"}\n)`,
						'data = data.merge(background, on="date", how="left")'
					].join('\n')
		});
	}

	list.push({
		id: 'look',
		title: 'Look at the data',
		text: 'Always plot the measurements before analysing them. Look for gaps, and for anything that seems out of place.',
		code: r ? `timePlot(data, pollutant = ${pollutant})` : `data.plot(x="date", y=${pollutant})`
	});

	if (isolating) {
		const column =
			settings.background === null ? null : control ? '"background"' : `"${settings.background}"`;
		const options = [
			...(column ? [`background${r ? ' = ' : '='}${dial('background', column)}`] : []),
			`deseason${r ? ' = ' : '='}${dial('deseason', bool[language](settings.deseason))}`,
			`deweather${r ? ' = ' : '='}${dial('deweather', bool[language](settings.deweather))}`
		];
		list.push({
			id: 'isolate',
			title: 'Isolate the signal',
			text:
				'Pollution rises and falls with the weather, the time of day and the time of year. ' +
				'This step fits a model of those patterns and takes them away, leaving the part of the signal that a change in emissions could explain.',
			code: r
				? `data$isolated <- isolateContribution(\n  data, ${pollutant},\n${options.map((o) => `  ${o}`).join(',\n')}\n)`
				: `data["isolated"] = aqeval.isolate_contribution(\n    data, ${pollutant},\n${options.map((o) => `    ${o},`).join('\n')}\n)`
		});
	}

	const averaging = dial('averaging', `"${settings.averaging}"`);
	list.push({
		id: 'average',
		title: 'Average the data',
		text:
			'Hourly values are noisy. Averaging them into longer blocks smooths out the short-lived ups and downs, ' +
			'and gives the next step fewer points to work through.',
		code: r
			? `data_avg <- timeAverage(data, avg.time = ${averaging})`
			: `data_avg = aqeval.time_average(data, ${averaging})`
	});

	const h = dial('h', String(settings.h));
	list.push({
		id: 'find',
		title: 'Find the break points',
		text:
			'A break point is a moment where the average level of the series shifts. ' +
			'h sets the shortest stretch allowed between two breaks, as a fraction of the whole series: a smaller h can find more breaks, closer together, and takes longer. ' +
			'0.3 is a quick first look. For 8-hour data we would normally use 0.12.',
		code: r
			? `breaks <- findBreakPoints(data_avg, ${target}, h = ${h})\nbreaks`
			: `breaks = aqeval.find_break_points(\n    data_avg, ${target}, h=${h}\n)\nbreaks`
	});

	list.push({
		id: 'quantify',
		title: 'Measure the changes',
		text:
			'Real changes rarely happen in a single day. This step fits a line through each stretch of the series and reports when each change started, when it finished, and how big it was.',
		code: r
			? `result <- quantBreakSegments(data_avg, ${target}, breaks)\nresult$report`
			: `result = aqeval.quant_break_segments(\n    data_avg, ${target}, breaks=breaks\n)\nresult["report"]`
	});

	return list;
}

/**
 * The whole tutorial as one script, with a comment heading each step.
 * @param {Settings} settings
 * @param {Language} language
 */
export function script(settings, language) {
	return steps(settings, language)
		.map((step, index) => `# ${index + 1}. ${step.title}\n${step.code}`)
		.join('\n\n');
}
