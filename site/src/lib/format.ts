/** Formatting for dates, numbers and the names of things. */

const DAY_MS = 86_400_000;

/** Timestamps from the pipeline are UTC wall-clock times without a zone. */
export function toMs(iso: string): number {
	return Date.parse(iso.length <= 10 ? `${iso}T00:00Z` : `${iso}Z`);
}

const dateFormat = new Intl.DateTimeFormat('en-GB', {
	day: 'numeric',
	month: 'short',
	year: 'numeric',
	timeZone: 'UTC'
});

const monthFormat = new Intl.DateTimeFormat('en-GB', { month: 'short', year: 'numeric', timeZone: 'UTC' });

export const formatDate = (value: string | number) =>
	dateFormat.format(typeof value === 'number' ? value : toMs(value));

export const formatMonth = (value: string | number) =>
	monthFormat.format(typeof value === 'number' ? value : toMs(value));

/** Whole days between two timestamps; positive when `to` is later. */
export const daysBetween = (from: string | number, to: string | number) =>
	Math.round(
		((typeof to === 'number' ? to : toMs(to)) - (typeof from === 'number' ? from : toMs(from))) / DAY_MS
	);

/** A length of time in the unit a person would use for it. */
export function formatDuration(days: number): string {
	const whole = Math.abs(days);
	if (whole < 1) return `${Math.round(whole * 24)} hours`;
	if (Math.round(whole) === 1) return '1 day';
	if (whole < 60) return `${Math.round(whole)} days`;
	if (whole < 730) return `${Math.round(whole / 30.44)} months`;
	return `${(whole / 365.25).toFixed(1)} years`;
}

export function formatNumber(value: number | null | undefined, digits = 1): string {
	if (value === null || value === undefined || Number.isNaN(value)) return '–';
	return value.toFixed(digits);
}

/** A change, with an explicit sign so a rise is never mistaken for a fall. */
export function formatSigned(value: number | null | undefined, digits = 1, suffix = ''): string {
	if (value === null || value === undefined || Number.isNaN(value)) return '–';
	// Round first, so that a value that rounds to nothing is not given a sign.
	const rounded = Number(value.toFixed(digits));
	const sign = rounded > 0 ? '+' : rounded < 0 ? '−' : '';
	return `${sign}${Math.abs(rounded).toFixed(digits)}${suffix}`;
}

export const UNIT = 'µg/m³';

export const POLLUTANT_LABEL: Record<string, string> = { no2: 'NO₂', nox: 'NOx' };

export const POLLUTANT_NAME: Record<string, string> = {
	no2: 'nitrogen dioxide',
	nox: 'nitrogen oxides'
};

export const AVERAGING_LABEL: Record<string, string> = {
	'8 hour': '8 hours',
	day: '1 day',
	'7 day': '1 week'
};

/** One line naming a set-up, for the comparison table. */
export function describeIsolation(settings: {
	deseason: boolean;
	deweather: boolean;
	background: string | null;
}): string {
	const parts = [
		settings.deseason ? 'cycles removed' : null,
		settings.deweather ? 'wind removed' : null,
		settings.background === 'air_temp'
			? 'air temperature'
			: settings.background
				? `background ${settings.background}`
				: null
	].filter(Boolean);
	return parts.length ? parts.join(', ') : 'as measured';
}
