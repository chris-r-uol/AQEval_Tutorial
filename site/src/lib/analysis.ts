/**
 * Turning a precomputed result into the things the page says about it.
 *
 * Nothing here re-does the analysis. It reads AQEval's break points and
 * segment report and works out the plain statements a reader wants: how many
 * changes, how big, and how close to the events they are being compared with.
 */
import { daysBetween, toMs } from './format';
import type { EventMarker, Run, Segment, Series } from './types';

/** A change smaller than this is the trend drifting, not a finding. */
export const MATERIAL_PERCENT = 5;

/** A change completed within this many days reads as a step rather than a slope. */
export const ABRUPT_DAYS = 60;

/** A series as [time, value] pairs, without the gaps. */
export function points(series: Series): Array<[number, number]> {
	const start = toMs(series.start);
	const step = series.stepHours * 3_600_000;
	const out: Array<[number, number]> = [];
	series.values.forEach((value, index) => {
		if (value !== null) out.push([start + index * step, value]);
	});
	return out;
}

export interface SegmentRow extends Segment {
	kind: 'fall' | 'rise' | 'steady';
	pace: 'abrupt' | 'gradual';
	material: boolean;
}

/** The segment report, with each row described in words. */
export function describeSegments(run: Run): SegmentRow[] {
	return run.segments.map((segment) => {
		const material = segment.percent !== null && Math.abs(segment.percent) >= MATERIAL_PERCENT;
		return {
			...segment,
			material,
			kind: !material ? 'steady' : (segment.change ?? 0) < 0 ? 'fall' : 'rise',
			pace: (segment.days ?? Infinity) <= ABRUPT_DAYS ? 'abrupt' : 'gradual'
		};
	});
}

/** Events are only tied to a change that happened within this many days. */
export const NEARBY_DAYS = 183;

export interface EventReading {
	event: EventMarker;
	/** How the event sits against the changes AQEval found. */
	relation: 'during' | 'near' | 'steady' | 'none';
	/** The change the event fell inside, or the nearest one. */
	row: SegmentRow | null;
	/** For 'near': days from the event to the change; negative when the change
	 *  finished before the event. */
	days: number | null;
}

/**
 * Where each event sits against the changes in the report.
 *
 * A date can fall inside a change, close to one, or in a stretch where the
 * trend was flat. Which of those it is — and not just whether a break point
 * landed near it — is what the reader is trying to find out.
 */
export function readEvents(rows: SegmentRow[], events: EventMarker[]): EventReading[] {
	const changes = rows.filter((row) => row.material);
	return events.map((event) => {
		const at = toMs(event.date);
		const during = changes.find((row) => toMs(row.from) <= at && at <= toMs(row.to));
		if (during) return { event, relation: 'during', row: during, days: 0 };

		let nearest: SegmentRow | null = null;
		let gap = Infinity;
		for (const row of changes) {
			// Days from the event to the nearer end of the change.
			const days = at < toMs(row.from) ? daysBetween(at, row.from) : -daysBetween(row.to, at);
			if (Math.abs(days) < Math.abs(gap)) {
				nearest = row;
				gap = days;
			}
		}
		if (nearest && Math.abs(gap) <= NEARBY_DAYS) return { event, relation: 'near', row: nearest, days: gap };

		const around = rows.find((row) => toMs(row.from) <= at && at <= toMs(row.to));
		return { event, relation: around ? 'steady' : 'none', row: around ?? null, days: null };
	});
}

export interface Summary {
	/** How many segments carry a change big enough to be a finding. */
	changes: number;
	/** The segment with the largest percentage change. */
	largest: SegmentRow | null;
	/** Fitted concentration at the start and end of the whole series. */
	start: number | null;
	end: number | null;
	overallPercent: number | null;
}

export function summarise(run: Run): Summary {
	const rows = describeSegments(run);
	const largest = rows
		.filter((row) => row.percent !== null)
		.reduce<SegmentRow | null>(
			(best, row) => (!best || Math.abs(row.percent!) > Math.abs(best.percent!) ? row : best),
			null
		);
	const start = run.trend[0]?.[1] ?? null;
	const end = run.trend[run.trend.length - 1]?.[1] ?? null;
	return {
		changes: rows.filter((row) => row.material).length,
		largest,
		start,
		end,
		overallPercent: start !== null && end !== null && start !== 0 ? ((end - start) / start) * 100 : null
	};
}
