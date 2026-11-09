<script lang="ts">
	/**
	 * What the step just added did, and what to look for.
	 *
	 * Shown in place of the findings until the last step, when there are
	 * findings to show. Each stage gets a few plain facts read from the data on
	 * screen, and a few things to try.
	 */
	import { profile, spread } from '$lib/analysis';
	import type { Settings, Stage } from '$lib/code';
	import { AVERAGE_NOUN, AVERAGING_LABEL, UNIT, daysBetween, formatDate, formatDuration, formatHour, formatNumber, toMs } from '$lib/format';
	import type { EventMarker, Run, Series } from '$lib/types';

	interface Props {
		stage: Stage;
		settings: Settings;
		/** The hourly measurements, on the first stage. */
		hourly: Series | null;
		/** The series the analysis runs on, averaged: isolated, or as measured. */
		averaged: Series | null;
		/** The measurements averaged the same way, to compare with. */
		measured: Series | null;
		isolating: boolean;
		run: Run | null;
		events: EventMarker[];
		pollutant: string;
		siteName: (code: string) => string;
		years: [number, number];
		language: 'r' | 'python';
	}

	let { stage, settings, hourly, averaged, measured, isolating, run, events, pollutant, siteName, years, language }: Props =
		$props();

	const count = (n: number) => n.toLocaleString('en-GB');

	let facts = $derived(hourly ? profile(hourly) : null);
	let span = $derived(years[1] - years[0] + 1);
	let caz = $derived(events.find((event) => event.id === 'caz'));

	let valid = $derived(averaged ? averaged.values.filter((value) => value !== null).length : 0);
	let blockHours = $derived(averaged?.stepHours ?? 0);
	let block = $derived(AVERAGE_NOUN[settings.averaging] ?? settings.averaging);

	// What the isolation step takes out, in words.
	let removed = $derived(
		[
			settings.deseason ? 'the daily and yearly cycle' : null,
			settings.deweather ? 'the wind' : null,
			settings.background === 'air_temp'
				? 'air temperature'
				: settings.background
					? `${pollutant} at ${siteName(settings.background)}`
					: null
		].filter((part): part is string => part !== null)
	);
	const list = (parts: string[]) =>
		parts.length < 2 ? parts.join('') : `${parts.slice(0, -1).join(', ')} and ${parts[parts.length - 1]}`;

	let before = $derived(measured ? spread(measured) : null);
	let after = $derived(averaged ? spread(averaged) : null);
	let smaller = $derived(before && after ? Math.round((1 - after / before) * 100) : null);

	// The break point closest to the Clean Air Zone, if there is one.
	let nearest = $derived.by(() => {
		if (!run || !caz || !run.breaks.length) return null;
		const ranked = [...run.breaks].sort(
			(a, b) => Math.abs(daysBetween(caz.date, a.date)) - Math.abs(daysBetween(caz.date, b.date))
		);
		const point = ranked[0];
		const at = toMs(caz.date);
		return { point, days: daysBetween(caz.date, point.date), covers: toMs(point.lower) <= at && at <= toMs(point.upper) };
	});
</script>

<div class="notes">
	{#if stage === 'data'}
		{#if facts}
			<dl class="tiles">
				<div class="tile">
					<dt>Hours in the record</dt>
					<dd>{count(facts.count)}</dd>
					<dd class="tile__note">{span} years, one value an hour</dd>
				</div>
				<div class="tile">
					<dt>Hours with a measurement</dt>
					<dd>{Math.round((facts.valid / facts.count) * 100)}%</dd>
					<dd class="tile__note">the rest are gaps in the line</dd>
				</div>
				<div class="tile">
					<dt>Most hours lie between</dt>
					<dd>{formatNumber(facts.typical[0], 0)}–{formatNumber(facts.typical[1], 0)}</dd>
					<dd class="tile__note">{UNIT}, the middle 90% of hours</dd>
				</div>
				<div class="tile">
					<dt>Highest hour</dt>
					<dd>{formatNumber(facts.highest.value, 0)}</dd>
					<dd class="tile__note">{UNIT}, {formatHour(facts.highest.time)}</dd>
				</div>
			</dl>
			<p>
				This is {pollutant} as the monitor recorded it, with nothing done to it. One value every hour for
				{span} years is {count(facts.count)} values, and from one day to the next they swing between about
				{formatNumber(facts.typical[0], 0)} and {formatNumber(facts.typical[1], 0)}
				{UNIT}.
			</p>
			<p>
				{#if caz}The Clean Air Zone started on {formatDate(caz.date)}.{/if} Can you tell from this figure
				whether the level changed then? The changes the analysis looks for are a few {UNIT} over months or
				years, much smaller than the swings from hour to hour. The steps that follow are how to find them.
			</p>
		{:else}
			<p class="muted">No measurements to describe.</p>
		{/if}
		<h3>Things to try</h3>
		<ul class="try">
			<li>Drag the ends of the bar under the figure to zoom in to a single week. When in the day is it highest, and why?</li>
			<li>Zoom back out to the whole record. Which time of year is worst?</li>
			<li>
				Choose another site. How does the rural background at Ladybower compare with the roadside at
				Bradford Mayo Avenue?
			</li>
		</ul>
	{:else if stage === 'average'}
		{#if averaged}
			<dl class="tiles">
				<div class="tile">
					<dt>{block[0].toUpperCase() + block.slice(1)} averages</dt>
					<dd>{count(valid)}</dd>
					<dd class="tile__note">from {count(averaged.values.length * blockHours)} hours</dd>
				</div>
				<div class="tile">
					<dt>Each point stands for</dt>
					<dd>{AVERAGING_LABEL[settings.averaging]}</dd>
					<dd class="tile__note">the mean of the hours measured in it</dd>
				</div>
			</dl>
		{/if}
		<p>
			Each point is the average of the measurements in one {block === 'daily' ? 'day' : block === 'weekly' ? 'week' : `${block} block`}. Averaging smooths out the short
			spikes, so slower changes are easier to see, and it leaves fewer values for the analysis to work
			through.
		</p>
		<p>
			The level still rises every winter and falls every summer, and it is higher in still weather than in
			wind. None of that is a change in what is emitted. The next step takes it out.
		</p>
		<h3>Things to try</h3>
		<ul class="try">
			<li>Switch between 8 hours, 1 day and 1 week in step 3. What can you see with a week that you could not with 8 hours? What do you lose?</li>
			<li>Go back to <em>Measurements</em> in the bar above the figure, and compare this with the hourly line.</li>
		</ul>
	{:else if stage === 'isolate'}
		{#if isolating}
			{#if smaller !== null && before !== null && after !== null}
				<dl class="tiles">
					<div class="tile">
						<dt>Spread of the measurements</dt>
						<dd>{formatNumber(before, 1)}</dd>
						<dd class="tile__note">{UNIT}, standard deviation</dd>
					</div>
					<div class="tile">
						<dt>Spread once isolated</dt>
						<dd>{formatNumber(after, 1)}</dd>
						<dd class="tile__note">
							{UNIT}, {smaller > 0 ? `${smaller}% less` : smaller < 0 ? `${-smaller}% more` : 'no change'}
						</dd>
					</div>
				</dl>
			{/if}
			<p>
				AQEval fitted a model of how {pollutant} at this site rises and falls with {list(removed)}, and took
				that away. What is left, in dark grey, is the isolated signal: the part of the measurements that the
				weather and the time of year do not explain, and so the part a change in emissions could have caused.
				The measurements are underneath in lighter grey.
			</p>
		{:else}
			<div class="alert">
				<p>
					<strong>Nothing is being removed.</strong> Tick at least one box in step 4, or choose something to
					account for, to isolate the signal.
				</p>
			</div>
		{/if}
		<p class="muted">
			In the code this step comes before the averaging: the model is fitted to the hourly measurements, and
			the isolated signal is then averaged. From here on, R and Python give slightly different answers,
			because they fit the model with different libraries. The page is showing the
			{language === 'r' ? 'R' : 'Python'} answers.
		</p>
		<h3>Things to try</h3>
		<ul class="try">
			<li>Untick one box at a time in step 4. Which makes the most difference to the spread?</li>
			<li>
				Choose a background site under <em>Also account for</em>. That takes out whatever the two sites have
				in common, such as changes across the whole region.
			</li>
		</ul>
	{:else if stage === 'breaks'}
		{#if run && run.breaks.length}
			<p>
				AQEval treats the series as a run of flat stretches and asks where one stretch ends and the next
				begins. With <code>h = {settings.h}</code> it found
				<strong>{run.breaks.length} break point{run.breaks.length === 1 ? '' : 's'}</strong>. Each is a
				dashed line on the figure, inside a shaded band: the range of dates it could fall in.
			</p>
			<ul class="points">
				{#each run.breaks as point (point.date)}
					<li>
						<strong>{formatDate(point.date)}</strong>
						<span class="muted">somewhere between {formatDate(point.lower)} and {formatDate(point.upper)}</span>
					</li>
				{/each}
			</ul>
			{#if nearest && caz}
				<p>
					{#if nearest.covers}
						The Clean Air Zone's start, {formatDate(caz.date)}, is inside the range of the break point on
						{formatDate(nearest.point.date)}.
					{:else}
						The break point nearest the Clean Air Zone is on {formatDate(nearest.point.date)},
						{formatDuration(Math.abs(nearest.days))}
						{nearest.days > 0 ? 'after' : 'before'} it started.
					{/if}
				</p>
			{/if}
		{:else if run}
			<div class="alert">
				<p>
					<strong>No break points at <code>h = {settings.h}</code>.</strong> With breaks forced this far
					apart, the series is best described as one flat stretch. Try a smaller <code>h</code> in step 5.
				</p>
			</div>
		{/if}
		<p>
			A break point says roughly when the level shifted. It does not say by how much, or whether the shift
			took a day or a year. That is the next step.
		</p>
		<h3>Things to try</h3>
		<ul class="try">
			<li>
				Make the search more sensitive: move <code>h</code> in step 5 down from 0.3 towards 0.12. Which break
				points stay, and which are new?
			</li>
			<li>Turn the isolation off again in step 4. Do the break points move?</li>
		</ul>
	{/if}
</div>

<style>
	.notes {
		font-size: 0.9375rem;
	}

	.notes p {
		max-width: 46rem;
	}

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

	.try,
	.points {
		margin: 0;
		padding-left: var(--space-5);
	}

	.try li,
	.points li {
		margin-bottom: var(--space-1);
	}
</style>
