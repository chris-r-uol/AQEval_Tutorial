/** Shapes of the JSON written by `pipeline/build_data.py`. */
import type { Language, Settings } from './code';

export interface PollutantStats {
	/** Annual mean concentration by year. */
	annual: Record<string, number | null>;
	/** Share of hours with a valid measurement, by year, in percent. */
	capture: Record<string, number | null>;
	/** Mean over the twelve months before the Clean Air Zone started. */
	before: number | null;
	/** Mean over the twelve months after. */
	after: number | null;
	mean: number | null;
}

export interface Site {
	code: string;
	name: string;
	/** AURN classification, e.g. "Urban Traffic". */
	type: string;
	authority: string;
	lat: number;
	lon: number;
	insideCaz: boolean;
	kmFromBradford: number;
	/** The isolation settings precomputed for this site. */
	isolations: string[];
	stats: Record<string, PollutantStats>;
}

export interface EventMarker {
	id: string;
	date: string;
	label: string;
	/** A few words, for the label on the figure. */
	short: string;
}

export interface Manifest {
	generated: string;
	/** The software that produced each language's results. */
	engines: Record<Language, string>;
	years: [number, number];
	sites: Site[];
	pollutants: string[];
	controlSites: string[];
	averaging: string[];
	h: number[];
	maxBreaksToQuantify: number;
	events: EventMarker[];
	centre: [number, number];
}

/** A regular time series: a start, a step, and the values in order. */
export interface Series {
	start: string;
	stepHours: number;
	values: Array<number | null>;
}

export interface BreakPoint {
	date: string;
	/** Confidence interval on the break's date. */
	lower: string;
	upper: string;
}

/** One row of AQEval's break-segment report. */
export interface Segment {
	from: string;
	to: string;
	days: number | null;
	/** Fitted concentration at the start and end of the segment. */
	c0: number | null;
	c1: number | null;
	change: number | null;
	percent: number | null;
	/** Values AQEval's report left blank because the data has a gap at that
	 *  moment, read from the fitted line instead. */
	filled?: Array<'c0' | 'c1'>;
}

export type RunStatus = 'ok' | 'no_breaks' | 'too_many_breaks' | 'error';

export interface Run {
	status: RunStatus;
	message?: string;
	breaks: BreakPoint[];
	/** The fitted trend as [time in ms, value, standard error] at each bend.
	 *  The standard error is null where the fit could not give one. */
	trend: Array<[number, number, number | null]>;
	segments: Segment[];
}

export interface RunFile {
	site: string;
	pollutant: string;
	deseason: boolean;
	deweather: boolean;
	background: string | null;
	/** The model AQEval fitted to isolate the signal, or null for none. */
	formula: string | null;
	averaging: Record<
		string,
		{
			start: string;
			stepHours: number;
			/** The series the analysis ran on: isolated, or as measured. */
			values: Array<number | null>;
			runs: Record<string, Run>;
		}
	>;
}

/** What the map labels each site with. */
export type MapMetric = 'before' | 'after' | 'change';

/** A result kept for comparison while the dials move on. */
export interface Pin {
	id: number;
	label: string;
	colour: string;
	settings: Settings;
	/** The language whose results these are. */
	language: Language;
	run: Run;
}
