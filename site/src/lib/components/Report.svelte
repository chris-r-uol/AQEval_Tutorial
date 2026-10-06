<script lang="ts">
	/**
	 * What the analysis found, in words and numbers.
	 *
	 * The figure shows the shape; this says what it means. It leads with the
	 * headline numbers, then reads each event against the changes, then gives
	 * AQEval's own report table — the same rows `result$report` prints. The
	 * break points come last: they are how the analysis got there, not what it
	 * found.
	 */
	import { describeSegments, readEvents, summarise, MATERIAL_PERCENT, type SegmentRow } from '$lib/analysis';
	import { UNIT, formatDate, formatDuration, formatNumber, formatSigned } from '$lib/format';
	import type { EventMarker, Run } from '$lib/types';

	interface Props {
		run: Run;
		events: EventMarker[];
		pollutant: string;
		maxBreaks: number;
		/** The language shown in the code panel, for the notes that depend on it. */
		language: 'r' | 'python';
	}

	let { run, events, pollutant, maxBreaks, language }: Props = $props();

	let rows = $derived(describeSegments(run));
	let summary = $derived(summarise(run));
	let readings = $derived(readEvents(rows, events));

	const what = (row: SegmentRow) =>
		row.kind === 'steady' ? 'Steady' : `${row.pace === 'abrupt' ? 'Abrupt' : 'Gradual'} ${row.kind}`;

	const size = (row: SegmentRow) => `${formatNumber(Math.abs(row.percent ?? 0), 0)}%`;
</script>

<div class="report">
	{#if run.status === 'error'}
		<div class="alert alert--error">
			<p><strong>AQEval could not fit a trend with this set-up.</strong></p>
			<p>
				This usually means a break sits too close to the start or end of the series. Try a different
				sensitivity or averaging period.
			</p>
		</div>
	{:else}
		<dl class="tiles">
			<div class="tile">
				<dt>Changes of {MATERIAL_PERCENT}% or more</dt>
				<dd>{summary.changes}</dd>
			</div>
			<div class="tile">
				<dt>Largest change</dt>
				{#if summary.largest && summary.largest.material}
					<dd class:fall={summary.largest.kind === 'fall'} class:rise={summary.largest.kind === 'rise'}>
						{formatSigned(summary.largest.percent, 0, '%')}
					</dd>
					<dd class="tile__note">
						over {formatDuration(summary.largest.days ?? 0)}, from {formatDate(summary.largest.from)}
					</dd>
				{:else}
					<dd>–</dd>
					<dd class="tile__note">nothing of {MATERIAL_PERCENT}% or more</dd>
				{/if}
			</div>
			<div class="tile">
				<dt>Start to end of the series</dt>
				<dd>{formatSigned(summary.overallPercent, 0, '%')}</dd>
				<dd class="tile__note">
					{formatNumber(summary.start)} → {formatNumber(summary.end)}
					{UNIT}
				</dd>
			</div>
		</dl>

		{#if run.status === 'no_breaks'}
			<div class="alert">
				<p>
					<strong>No break points at this sensitivity.</strong> With breaks forced this far apart, the series
					is best described as one straight line. Try a smaller <code>h</code>.
				</p>
			</div>
		{:else if run.status === 'too_many_breaks'}
			<div class="alert">
				<p>
					<strong>{run.breaks.length} break points is too many to measure here.</strong> To measure the
					changes, AQEval tries nine starting points for every break, in every combination. For
					{run.breaks.length} breaks that is {(9 ** run.breaks.length).toLocaleString('en-GB')} model fits,
					which takes too long to have ready for every set-up. The break points are shown on the figure.
					Choose a larger <code>h</code> or a longer averaging period to measure the changes, or run this
					set-up yourself and wait.
				</p>
			</div>
		{/if}

		{#if rows.length}
			<h3>What happened around each event</h3>
			<ul class="events">
				{#each readings as reading (reading.event.id)}
					<li>
						<strong>{reading.event.label}</strong>
						<span class="muted">{formatDate(reading.event.date)}</span>
						<p>
							{#if reading.relation === 'during' && reading.row}
								This date is inside {reading.row.pace === 'abrupt' ? 'an abrupt' : 'a gradual'}
								{reading.row.kind} of {size(reading.row)}, which ran from {formatDate(reading.row.from)} to
								{formatDate(reading.row.to)}.
							{:else if reading.relation === 'near' && reading.row}
								The nearest change is {reading.row.pace === 'abrupt' ? 'an abrupt' : 'a gradual'}
								{reading.row.kind} of {size(reading.row)}
								{#if (reading.days ?? 0) > 0}
									that started {formatDuration(reading.days ?? 0)} later, on {formatDate(reading.row.from)}.
								{:else}
									that finished {formatDuration(reading.days ?? 0)} earlier, on {formatDate(reading.row.to)}.
								{/if}
							{:else if reading.relation === 'steady'}
								The trend was steady across this date: no change of {MATERIAL_PERCENT}% or more within six
								months.
							{:else}
								No fitted trend at this date.
							{/if}
						</p>
					</li>
				{/each}
			</ul>

			<h3>The report</h3>
			<p class="muted">
				One row for each stretch of the fitted trend. These are the rows that
				<code>{language === 'r' ? 'result$report' : 'result["report"]'}</code> prints.
			</p>
			<div class="table-scroll">
				<table>
					<thead>
						<tr>
							<th scope="col">Period</th>
							<th scope="col" class="num">Start → end</th>
							<th scope="col" class="num">Change</th>
							<th scope="col">What happened</th>
						</tr>
					</thead>
					<tbody>
						{#each rows as row (row.from)}
							<tr class:quiet={!row.material}>
								<td>
									{formatDate(row.from)} to {formatDate(row.to)}
									<span class="sub">{row.days === null ? '' : formatDuration(row.days)}</span>
								</td>
								<td class="num">
									{formatNumber(row.c0)}{row.filled?.includes('c0') ? '*' : ''} →
									{formatNumber(row.c1)}{row.filled?.includes('c1') ? '*' : ''}
								</td>
								<td class="num">
									{formatSigned(row.change)}
									<span class="sub">{formatSigned(row.percent, 0, '%')}</span>
								</td>
								<td>
									<span class="what" class:what--fall={row.kind === 'fall'} class:what--rise={row.kind === 'rise'}>
										{what(row)}
									</span>
								</td>
							</tr>
						{/each}
					</tbody>
				</table>
			</div>
			<p class="muted footnote">
				Concentrations are in {UNIT} of {pollutant}, read from the fitted trend. A change under
				{MATERIAL_PERCENT}% is shown as steady; one completed within 60 days as abrupt.
				{#if rows.some((row) => row.filled)}
					* The measurements have a gap at this moment, so AQEval's own report leaves the value blank, and
					the change with it. It is read here from the fitted line.
				{/if}
			</p>
		{/if}

		{#if run.breaks.length}
			<details class="breaks" open={rows.length === 0}>
				<summary>The {run.breaks.length} break point{run.breaks.length === 1 ? '' : 's'} from the step before</summary>
				<p class="muted">
					The search for break points asks a simpler question: where does the average level shift? Its
					answers are the starting guesses for the report above, not the result. The report then works out
					where each change really begins and ends, which is why the two rarely line up. Switch on
					<em>Break points</em> in the key of the figure to see them.
				</p>
				<ul>
					{#each run.breaks as point (point.date)}
						<li>
							<strong>{formatDate(point.date)}</strong>
							<span class="muted">between {formatDate(point.lower)} and {formatDate(point.upper)}</span>
						</li>
					{/each}
				</ul>
			</details>
		{/if}
	{/if}
</div>

<style>
	h3 {
		font-size: 0.9375rem;
		margin: var(--space-4) 0 var(--space-1);
	}

	.tiles {
		display: grid;
		grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
		gap: var(--space-3);
		margin: 0 0 var(--space-3);
	}

	.tile {
		padding: var(--space-3);
		background: var(--surface-sunken);
		border-radius: var(--radius);
	}

	.tile dt {
		font-size: 0.8125rem;
		color: var(--text-muted);
	}

	.tile dd {
		margin: 0;
		font-size: 1.75rem;
		font-weight: 700;
		line-height: 1.2;
		font-variant-numeric: tabular-nums;
	}

	.tile dd.tile__note {
		font-size: 0.8125rem;
		font-weight: 400;
		color: var(--text-muted);
		line-height: 1.35;
	}

	.fall {
		color: var(--colour-fall);
	}

	.rise {
		color: var(--colour-rise);
	}

	.events {
		list-style: none;
		margin: 0;
		padding: 0;
	}

	.breaks {
		margin-top: var(--space-4);
		font-size: 0.9375rem;
	}

	.breaks summary {
		color: var(--colour-brand);
		cursor: pointer;
		text-decoration: underline;
		width: fit-content;
	}

	.breaks p {
		margin: var(--space-2) 0;
	}

	.breaks ul {
		margin: 0;
		padding-left: var(--space-5);
	}

	.events li {
		padding: var(--space-2) 0 var(--space-2) var(--space-3);
		border-left: 3px solid var(--text);
		margin-bottom: var(--space-2);
		font-size: 0.9375rem;
	}

	.events p {
		margin: 2px 0 0;
	}

	.table-scroll {
		overflow-x: auto;
	}

	table {
		width: 100%;
		border-collapse: collapse;
		font-size: 0.875rem;
		font-variant-numeric: tabular-nums;
	}

	th,
	td {
		padding: var(--space-2);
		text-align: left;
		border-bottom: 1px solid var(--border-subtle);
		vertical-align: top;
	}

	th,
	.num,
	.what {
		white-space: nowrap;
	}

	/* The second line of a cell: the length of a period, or a change as a percentage. */
	.sub {
		display: block;
		color: var(--text-muted);
		font-size: 0.8125rem;
	}

	th {
		font-weight: 700;
		border-bottom: 2px solid var(--text);
	}

	.num {
		text-align: right;
	}

	tr.quiet td {
		color: var(--text-muted);
	}

	.what {
		display: inline-block;
		padding: 1px var(--space-2);
		border-radius: 999px;
		background: var(--surface-sunken);
		font-size: 0.8125rem;
		font-weight: 600;
	}

	.what--fall {
		background: var(--colour-fall-tint);
		color: var(--colour-fall);
	}

	.what--rise {
		background: var(--colour-rise-tint);
		color: var(--colour-rise);
	}

	.footnote {
		margin: var(--space-2) 0 0;
	}
</style>
