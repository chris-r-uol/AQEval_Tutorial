/**
 * Chart option builders.
 *
 * Kept out of the components so every chart shares one axis treatment, one
 * tooltip style and one palette.
 */
import type { EChartsOption, SeriesOption } from 'echarts';
import { formatDate, toMs } from './format';
import type { SegmentRow } from './analysis';
import type { EventMarker, Run, Site } from './types';

/** AQEval's own figures use grey data, a red trend and blue for the changes.
 *  Anyone who runs the code will see those colours, so the page uses them too. */
export const COLOURS = {
	measured: '#cfd2d3',
	analysed: '#8d9294',
	trend: '#d4351c',
	break: '#1d52b8',
	breakPoint: '#505a5f',
	event: '#0b0c0c'
};

const CHANGE_SERIES = 'Change';
const BREAK_POINT_SERIES = 'Break points';

/** Colours for pinned results, chosen to stay apart from the red trend. */
export const PIN_COLOURS = ['#5f2c6d', '#00703c', '#b58840', '#28a197'];

export const SITE_COLOURS: Record<string, string> = {
	'Urban Traffic': '#912b88',
	'Urban Background': '#28a197'
};

const AXIS_STYLE = {
	axisLine: { lineStyle: { color: '#b1b4b6' } },
	axisLabel: { color: '#505a5f', fontSize: 11 },
	splitLine: { lineStyle: { color: '#e5e5e5' } }
};

const TOOLTIP = {
	backgroundColor: 'rgba(11,12,12,0.92)',
	borderWidth: 0,
	textStyle: { color: '#ffffff', fontSize: 12 }
};

const TEXT = { fontFamily: 'arial, helvetica, sans-serif' };

export interface PinnedTrend {
	label: string;
	colour: string;
	trend: Run['trend'];
}

export interface FigureInput {
	/** The series the analysis ran on, as [time, value]. */
	analysed: Array<[number, number]>;
	analysedName: string;
	/** The measured series, drawn underneath when the analysis ran on an
	 *  isolated one. */
	measured: Array<[number, number]> | null;
	run: Run;
	/** The run's segment report, described. */
	rows: SegmentRow[];
	events: EventMarker[];
	pins: PinnedTrend[];
	unit: string;
	range: [number, number];
}

/**
 * A vertical range that holds the middle 99% of the points.
 *
 * A handful of extreme 8-hour values would otherwise set the scale, and the
 * trend, which moves by a few units over years, would be a flat line in the
 * middle of it. The figure's caption says that the outliers are not drawn.
 */
export function focusRange(values: Array<[number, number]>): [number, number] {
	const sorted = values.map(([, value]) => value).sort((a, b) => a - b);
	if (!sorted.length) return [0, 1];
	const low = sorted[Math.floor(0.005 * (sorted.length - 1))];
	const high = sorted[Math.ceil(0.995 * (sorted.length - 1))];
	return [Math.floor(low / 5) * 5, Math.ceil(high / 5) * 5];
}

/** The main figure: data, fitted trend, break points and events. */
export function figureOption(input: FigureInput): EChartsOption {
	const { analysed, analysedName, measured, run, rows, events, pins, unit, range } = input;
	const [low, high] = focusRange(analysed);

	const scatter = (name: string, data: Array<[number, number]>, colour: string, z: number): SeriesOption => ({
		name,
		// Points, not a line: these are averages of separate periods, and
		// joining them implies a continuity the measurements do not have.
		// Progressive rather than `large`, which does not cope with a time axis.
		type: 'scatter',
		symbolSize: data.length > 3000 ? 2.2 : data.length > 600 ? 3 : 4.5,
		itemStyle: { color: colour },
		data,
		// Points beyond the vertical range are left out, not piled on its edge.
		clip: true,
		z,
		progressive: 2000,
		progressiveThreshold: 2000,
		tooltip: { show: false }
	});

	const series: SeriesOption[] = [];
	if (measured) series.push(scatter('Measured', measured, COLOURS.measured, 1));
	series.push(scatter(analysedName, analysed, COLOURS.analysed, 2));

	for (const pin of pins) {
		series.push({
			name: `Pinned ${pin.label}`,
			type: 'line',
			symbol: 'none',
			lineStyle: { width: 2, color: pin.colour, type: [6, 4] },
			itemStyle: { color: pin.colour },
			data: pin.trend.map(([time, value]) => [time, value]),
			z: 4
		});
	}

	series.push({
		name: 'Fitted trend',
		type: 'line',
		symbol: 'none',
		lineStyle: { width: 2.6, color: COLOURS.trend },
		itemStyle: { color: COLOURS.trend },
		data: run.trend.map(([time, value]) => [time, value]),
		z: 6
	});

	// Blue marks where the trend changes, which is the result of the analysis.
	// A change that was over within weeks is shaded as a band, because it
	// happened across that whole window rather than on one day; where one long
	// stretch gives way to another there is a line. The marks hang off a series
	// with no data of its own, so the legend can switch them on and off.
	const bands = rows
		.filter((row) => row.pace === 'abrupt')
		.map((row) => [{ xAxis: toMs(row.from) }, { xAxis: toMs(row.to) }]);
	const boundaries = rows
		.slice(1)
		.filter((row, index) => row.pace !== 'abrupt' && rows[index].pace !== 'abrupt')
		.map((row) => ({ xAxis: toMs(row.from) }));

	series.push({
		name: CHANGE_SERIES,
		type: 'line',
		data: [],
		itemStyle: { color: COLOURS.break },
		lineStyle: { color: COLOURS.break },
		markArea: {
			silent: true,
			itemStyle: { color: 'rgba(29, 82, 184, 0.16)' },
			data: bands as never
		},
		markLine: {
			silent: true,
			symbol: 'none',
			label: { show: false },
			lineStyle: { color: COLOURS.break, width: 1.6, type: 'dashed' },
			data: boundaries
		}
	});

	// The break points from the step before are the starting guesses handed to
	// the segment fit, not its answer, and they rarely sit where the trend
	// bends. They are there to be looked at, but start switched off.
	series.push({
		name: BREAK_POINT_SERIES,
		type: 'line',
		data: [],
		itemStyle: { color: COLOURS.breakPoint },
		lineStyle: { color: COLOURS.breakPoint },
		markArea: {
			silent: true,
			// How sure the search is about each date.
			itemStyle: { color: 'rgba(80, 90, 95, 0.12)' },
			data: run.breaks.map((b) => [{ xAxis: toMs(b.lower) }, { xAxis: toMs(b.upper) }]) as never
		},
		markLine: {
			silent: true,
			symbol: 'none',
			label: { show: false },
			lineStyle: { color: COLOURS.breakPoint, width: 1.2, type: 'dotted' },
			data: run.breaks.map((b) => ({ xAxis: toMs(b.date) }))
		}
	});

	series.push({
		name: 'Event',
		type: 'line',
		data: [],
		itemStyle: { color: COLOURS.event },
		lineStyle: { color: COLOURS.event },
		markLine: {
			silent: true,
			symbol: 'none',
			lineStyle: { color: COLOURS.event, width: 1.2, type: 'solid' },
			label: {
				formatter: (params) => String((params.data as { name?: string }).name ?? ''),
				// Upright, above the top of the line, clear of the data.
				position: 'end',
				color: COLOURS.event,
				fontSize: 11,
				fontWeight: 'bold'
			},
			data: events.map((event) => ({ xAxis: toMs(event.date), name: event.short }))
		}
	});

	return {
		textStyle: TEXT,
		animation: false,
		grid: { top: 58, right: 18, bottom: 58, left: 50 },
		legend: {
			// One row that pages sideways on a narrow screen, rather than
			// wrapping onto the top of the plot.
			type: 'scroll',
			top: 0,
			left: 0,
			itemHeight: 8,
			itemWidth: 16,
			textStyle: { fontSize: 12, color: '#505a5f' },
			// With no segments to show, the break points are all there is.
			selected: { [BREAK_POINT_SERIES]: rows.length === 0 }
		},
		tooltip: {
			...TOOLTIP,
			trigger: 'axis',
			axisPointer: { type: 'line', lineStyle: { color: '#505a5f' } },
			formatter: (params) => {
				const list = Array.isArray(params) ? params : [params];
				if (!list.length) return '';
				const time = (list[0].value as number[])[0];
				const rows = list
					.filter((item) => item.seriesType === 'line')
					.map((item) => `${item.marker}${item.seriesName}: <b>${(item.value as number[])[1].toFixed(1)}</b> ${unit}`);
				return [formatDate(time), ...rows].join('<br>');
			}
		},
		xAxis: { type: 'time', min: range[0], max: range[1], ...AXIS_STYLE, splitLine: { show: false } },
		yAxis: {
			type: 'value',
			name: unit,
			nameTextStyle: { color: '#505a5f', fontSize: 11, align: 'left' },
			min: low,
			max: high,
			...AXIS_STYLE
		},
		// A slider only. Zooming on the mouse wheel would capture the wheel as
		// the reader scrolls the page past the figure.
		dataZoom: [
			{ type: 'slider', height: 18, bottom: 8, filterMode: 'none', borderColor: '#e5e5e5', showDetail: false }
		],
		series
	};
}

/** Annual means for every site, so the sites can be compared at a glance. */
export function annualOption(
	sites: Site[],
	pollutant: string,
	selected: string,
	unit: string,
	events: EventMarker[]
): EChartsOption {
	const years = Object.keys(sites[0].stats[pollutant].annual);
	return {
		textStyle: TEXT,
		animation: false,
		grid: { top: 52, right: 16, bottom: 26, left: 44 },
		legend: {
			type: 'scroll',
			top: 0,
			left: 0,
			itemHeight: 8,
			itemWidth: 16,
			textStyle: { fontSize: 11, color: '#505a5f' }
		},
		tooltip: { ...TOOLTIP, trigger: 'axis', valueFormatter: (value) => `${Number(value).toFixed(1)} ${unit}` },
		xAxis: { type: 'category', data: years, ...AXIS_STYLE, splitLine: { show: false } },
		yAxis: {
			type: 'value',
			name: unit,
			nameTextStyle: { color: '#505a5f', fontSize: 11, align: 'left' },
			...AXIS_STYLE
		},
		series: sites.map((site, index) => ({
			name: site.name,
			type: 'line',
			symbolSize: site.code === selected ? 7 : 5,
			// The selected site is drawn heavier; line style separates the two
			// sites of each type so colour is not doing all the work.
			lineStyle: {
				width: site.code === selected ? 3.2 : 1.6,
				type: sites.findIndex((s) => s.type === site.type) === index ? 'solid' : 'dashed'
			},
			color: SITE_COLOURS[site.type] ?? '#505a5f',
			z: site.code === selected ? 5 : 2,
			data: years.map((year) => site.stats[pollutant].annual[year]),
			// The events are drawn once, on the first series.
			markLine:
				index === 0
					? {
							silent: true,
							symbol: 'none',
							label: { show: false },
							lineStyle: { color: '#0b0c0c', width: 1, type: 'dotted' },
							data: events.map((event) => ({ xAxis: event.date.slice(0, 4) }))
						}
					: undefined
		}))
	};
}
