export type Language = 'r' | 'python';

export interface Settings {
	site: string;
	pollutant: string;
	deseason: boolean;
	deweather: boolean;
	background: string | null;
	averaging: string;
	h: number;
}

export interface Step {
	id: string;
	title: string;
	text: string;
	code: string;
}

export const YEARS: [number, number];
export const DEFAULTS: Settings;
export function plain(code: string): string;
export function segments(code: string): Array<{ text: string; dial: keyof Settings | null }>;
export function isolates(settings: Settings): boolean;
export function usesControlSite(settings: Settings): boolean;
export function steps(settings: Settings, language: Language): Step[];
export function script(settings: Settings, language: Language): string;
