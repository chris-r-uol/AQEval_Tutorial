<script lang="ts">
	/**
	 * The experiment.
	 *
	 * Dials on the left, results in the middle, the code that produces them on
	 * the right. Every result was computed ahead of time by the real AQEval
	 * analysis (see `pipeline/build_data.py`); the page looks up the one that
	 * matches the dials.
	 *
	 * A reader starts with the measurements alone and adds the analysis a step
	 * at a time. The figure, the notes under it and the code show only the
	 * steps added so far.
	 */
	import { onMount, untrack } from 'svelte';
	import { afterNavigate, replaceState } from '$app/navigation';
	import type { Feature, Polygon } from 'geojson';
	import { describeSegments, points, timeline } from '$lib/analysis';
	import { annualOption, figureOption, PIN_COLOURS, ratioOption } from '$lib/charts';
	import { isolates, reached, type Language, type Settings, type Stage } from '$lib/code';
	import { CODESPACE_URL, REPOSITORY_URL } from '$lib/config';
	import { loadBoundary, loadHourly, loadManifest, loadMeasured, loadRuns } from '$lib/data';
	import { experiment } from '$lib/experiment.svelte';
	import {
		AVERAGE_NOUN,
		AVERAGING_LABEL,
		POLLUTANT_LABEL,
		UNIT,
		describeIsolation,
		formatDate
	} from '$lib/format';
	import { STAGE_ACTION, STAGE_TEASER } from '$lib/stages';
	import type { Manifest, MapMetric, Pin, RunFile, Series } from '$lib/types';
	import Chart from '$lib/components/Chart.svelte';
	import CodePanel from '$lib/components/CodePanel.svelte';
	import Compare from '$lib/components/Compare.svelte';
	import Controls from '$lib/components/Controls.svelte';
	import Report from '$lib/components/Report.svelte';
	import Segmented from '$lib/components/Segmented.svelte';
	import SiteMap from '$lib/components/SiteMap.svelte';
	import StageNotes from '$lib/components/StageNotes.svelte';
	import Stepper from '$lib/components/Stepper.svelte';
	import Tasks from '$lib/components/Tasks.svelte';

	let manifest = $state.raw<Manifest | null>(null);
	let boundary = $state.raw<Feature<Polygon> | null>(null);
	let failure = $state<string | null>(null);

	/** What is on screen: one set of dials and the data that goes with it.
	 *  Replaced whole, so the figure never mixes one set-up's settings with
	 *  another's data while a file is still on its way. */
	interface View {
		/** The settings in force at this stage: before the isolation step is
		 *  added, nothing is isolated, whatever its dials say. */
		settings: Settings;
		language: Language;
		stage: Stage;
		/** The hourly measurements, for the first stage. */
		hourly: Series | null;
		/** The results, once there is an analysis to show. */
		runs: RunFile | null;
		/** The same site and pollutant with nothing isolated. */
		measured: RunFile | null;
	}
	let view = $state.raw<View | null>(null);
	let loading = $state(false);

	const allowed = (m: Manifest) => ({
		sites: m.sites.map((site) => site.code),
		pollutants: m.pollutants,
		averaging: m.averaging,
		h: m.h,
		controls: m.controlSites
	});

	onMount(() => {
		loadManifest()
			.then((loaded) => {
				experiment.fromHash(location.hash, allowed(loaded));
				manifest = loaded;
			})
			.catch((error) => (failure = String(error.message ?? error)));
		// The boundary is decoration; the page works without it.
		loadBoundary()
			.then((loaded) => (boundary = loaded))
			.catch(() => {});

		const onHashChange = () => manifest && experiment.fromHash(location.hash, allowed(manifest));
		window.addEventListener('hashchange', onHashChange);
		return () => window.removeEventListener('hashchange', onHashChange);
	});

	// Fetch the results for the current dials and stage. Files already fetched
	// resolve at once, so only a new site, pollutant or isolation setting waits.
	$effect(() => {
		if (!manifest) return;
		const settings = experiment.effective;
		const language = experiment.language;
		const stage = experiment.stage;
		const analysing = reached(stage, 'average');
		let cancelled = false;
		loading = true;
		Promise.all([
			analysing ? null : loadHourly(settings),
			analysing ? loadRuns(settings, language) : null,
			analysing ? loadMeasured(settings, language) : null
		])
			.then(([hourly, runs, measured]) => {
				if (cancelled) return;
				view = { settings, language, stage, hourly, runs, measured };
				failure = null;
				loading = false;
			})
			.catch((error) => {
				if (cancelled) return;
				failure = String(error.message ?? error);
				loading = false;
			});
		// Fetch what the next step needs while this one is being read, so that
		// adding it is instant. A failure here is reported when it is needed.
		const ahead = untrack(() => experiment.settings);
		if (!analysing) loadMeasured(settings, language).catch(() => {});
		else if (!reached(stage, 'isolate')) loadRuns(ahead, language).catch(() => {});
		return () => {
			cancelled = true;
		};
	});

	// Keep the address in step with the dials, so the page can be linked to.
	// SvelteKit's router owns the address and is not ready until the first
	// navigation has finished.
	let routerReady = $state(false);
	afterNavigate(() => {
		routerReady = true;
	});

	$effect(() => {
		if (!manifest || !routerReady) return;
		const hash = experiment.toHash();
		const clean = experiment.isDefault && experiment.language === 'r';
		if ((clean && location.hash) || (!clean && location.hash !== hash)) {
			replaceState(clean ? location.pathname + location.search : hash, {});
		}
	});

	let siteName = (code: string) => manifest?.sites.find((site) => site.code === code)?.name ?? code;

	let averaged = $derived(view?.runs ? view.runs.averaging[view.settings.averaging] : null);
	let measuredAveraged = $derived(view?.measured ? view.measured.averaging[view.settings.averaging] : null);
	// The run holds the break points and the changes, so it only counts once
	// the reader has asked for them.
	let run = $derived(
		averaged && view && reached(view.stage, 'breaks') ? averaged.runs[String(view.settings.h)] : null
	);
	let measuring = $derived(view?.stage === 'changes');
	let isolating = $derived(view ? isolates(view.settings) : false);
	let pollutant = $derived(view ? POLLUTANT_LABEL[view.settings.pollutant] : '');

	let showMeasured = $state(true);

	// On a wide screen the code has a column of its own. Otherwise it follows
	// the figure, so the two are still read together. Rendered in one place or
	// the other rather than moved with CSS, because the reading order differs.
	const WIDE = '(min-width: 1360px)';
	let wide = $state(false);

	onMount(() => {
		const query = window.matchMedia(WIDE);
		wide = query.matches;
		const onChange = (event: MediaQueryListEvent) => (wide = event.matches);
		query.addEventListener('change', onChange);
		return () => query.removeEventListener('change', onChange);
	});

	// --- Pins ------------------------------------------------------------------

	const MAX_PINS = PIN_COLOURS.length;
	let pins = $state.raw<Pin[]>([]);
	let nextPinId = 1;

	const sameSettings = (a: Settings, b: Settings) =>
		(Object.keys(a) as Array<keyof Settings>).every((key) => a[key] === b[key]);

	/** Whether a pin is the result on screen: same dials, same language. */
	const onScreen = (p: Pin) => !!view && p.language === view.language && sameSettings(p.settings, view.settings);

	let alreadyPinned = $derived(pins.some(onScreen));

	function pin() {
		if (!view || !run || !measuring || pins.length >= MAX_PINS || alreadyPinned) return;
		// Reuse the first letter that is free, so labels stay A to D.
		const index = PIN_COLOURS.findIndex((_, i) => !pins.some((p) => p.label === 'ABCD'[i]));
		pins = [
			...pins,
			{
				id: nextPinId++,
				label: 'ABCD'[index],
				colour: PIN_COLOURS[index],
				settings: view.settings,
				language: view.language,
				run
			}
		];
	}

	// --- Figure ----------------------------------------------------------------

	let rows = $derived(run && measuring ? describeSegments(run) : []);

	let figure = $derived.by(() => {
		if (!view || !manifest) return null;
		const range: [number, number] = [Date.UTC(manifest.years[0], 0, 1), Date.UTC(manifest.years[1] + 1, 0, 1)];
		if (!reached(view.stage, 'average')) {
			if (!view.hourly) return null;
			return figureOption({
				analysed: points(view.hourly),
				analysedName: 'Hourly measurements',
				hourly: timeline(view.hourly),
				unit: UNIT,
				range
			});
		}
		if (!averaged || !measuredAveraged) return null;
		const measured = points(measuredAveraged);
		const isolated = isolating ? points(averaged) : null;
		return figureOption({
			analysed: isolated ?? measured,
			analysedName: isolated ? 'Isolated signal' : 'Measured',
			measured: isolated && showMeasured ? measured : null,
			breaks: run?.breaks ?? null,
			trend: measuring && run ? run.trend : null,
			rows,
			// The events come in with the break points, which is when the
			// question becomes when the level changed.
			events: run ? manifest.events : [],
			// A pin of the set-up on screen would only draw over the trend.
			pins: measuring
				? pins
						.filter((p) => !onScreen(p))
						.map((p) => ({ label: p.label, colour: p.colour, trend: p.run.trend }))
				: [],
			unit: UNIT,
			range
		});
	});

	let averagingLabel = $derived(view ? AVERAGING_LABEL[view.settings.averaging] : '');

	let figureDescription = $derived.by(() => {
		if (!view) return 'Loading';
		const where = `${pollutant} at ${siteName(view.settings.site)}`;
		switch (view.stage) {
			case 'data':
				return `Hourly ${where}, as measured, ${manifest?.years[0]} to ${manifest?.years[1]}.`;
			case 'average':
				return `${where}, averaged to ${averagingLabel}.`;
			case 'isolate':
				return `${where}, ${isolating ? 'with the weather and seasons removed, ' : ''}averaged to ${averagingLabel}.`;
			case 'breaks':
				return `${where}, averaged to ${averagingLabel}, with ${run?.breaks.length ?? 0} break points marked. The notes below the figure give their dates.`;
			default:
				return (
					`${where}, averaged to ${averagingLabel}, ` +
					`with the fitted trend in ${rows.length || 1} stretch${rows.length > 1 ? 'es' : ''}. ` +
					'The findings below the figure give the dates and sizes of the changes.'
				);
		}
	});

	const NOTES_TITLE: Record<Stage, string> = {
		data: 'What you are looking at',
		average: 'What averaging does',
		isolate: 'What isolating the signal does',
		breaks: 'Where the level shifts',
		changes: 'What the analysis found'
	};

	// --- Map and site comparison -------------------------------------------------

	let metric = $state<MapMetric>('change');
	let caz = $derived(manifest?.events.find((event) => event.id === 'caz'));

	let annual = $derived(
		manifest
			? annualOption(manifest.sites, experiment.settings.pollutant, experiment.settings.site, UNIT, manifest.events)
			: null
	);

	let ratio = $derived(manifest ? ratioOption(manifest.sites, experiment.settings.site, manifest.events) : null);

	let selectedSite = $derived(manifest?.sites.find((site) => site.code === experiment.settings.site));
	let captureYears = $derived(
		selectedSite ? Object.entries(selectedSite.stats[experiment.settings.pollutant].capture) : []
	);
</script>

<section class="intro">
	<h1>Did Bradford's Clean Air Zone change the air?</h1>
	<p class="intro__lede">
		Bradford began charging the most polluting commercial vehicles to enter the city on 26 September 2022.
		This page is the experiment an air quality analyst would run to find out what happened next: take the
		measurements, remove the weather and the seasons, and ask whether the level changed, when, and by how
		much. The method is <a href="https://karlropkins.github.io/AQEval/" target="_blank" rel="noopener">AQEval</a>.
	</p>
	<p class="intro__how">
		<strong>Start with the measurements</strong>, as the monitor recorded them. Then add the analysis one
		step at a time, and watch the figure, the notes and the code build up together. When you are ready to
		run the code for real, <a href={CODESPACE_URL} target="_blank" rel="noopener">open the workspace</a>.
	</p>
</section>

{#if failure && !manifest}
	<div class="alert alert--error load-error" role="alert">
		<p><strong>The results could not be loaded.</strong> {failure}</p>
		<p>Check your connection and reload the page.</p>
	</div>
{:else if !manifest}
	<div class="bench" aria-busy="true">
		<div class="card skeleton" style:height="520px"></div>
		<div class="card skeleton" style:height="520px"></div>
	</div>
{:else}
	<div class="bench" id="experiment">
		<aside class="bench__controls card" aria-label="Experiment set-up">
			<Controls {manifest} />
		</aside>

		<div class="bench__main">
			<section class="card" aria-labelledby="figure-title">
				<Stepper />
				<div class="card__head">
					<div>
						<h2 id="figure-title" class="card__title">
							{#if view}{pollutant} at {siteName(view.settings.site)}{:else if failure}No result{:else}Loading…{/if}
						</h2>
						{#if view}
							<p class="muted setup">
								{#if view.stage === 'data'}
									Hourly measurements, {manifest.years[0]} to {manifest.years[1]} · nothing averaged or removed
								{:else if view.stage === 'average'}
									As measured · averaged to {averagingLabel}
								{:else}
									{describeIsolation(view.settings)} · averaged to {averagingLabel}
									{#if run}· <code>h = {view.settings.h}</code>{/if}
									· computed in {view.language === 'r' ? 'R' : 'Python'}
								{/if}
							</p>
						{/if}
					</div>
					{#if measuring}
						<button
							class="btn btn--secondary btn--small"
							onclick={pin}
							disabled={!run || alreadyPinned || pins.length >= MAX_PINS}
							title={pins.length >= MAX_PINS ? `Remove a pin first: ${MAX_PINS} is the most the figure can show` : ''}
						>
							{alreadyPinned ? 'Pinned' : 'Pin this result'}
						</button>
					{/if}
				</div>

				<div class="figure" class:figure--loading={loading}>
					{#if figure}
						<Chart option={figure} height="400px" description={figureDescription} />
					{:else if !failure}
						<div class="skeleton" style:height="400px"></div>
					{/if}
				</div>

				{#if view && view.stage === 'data' && view.hourly}
					<p class="muted caption">
						The grey line is every hourly measurement, as the monitor recorded it. Where it stops, the
						monitor has no valid measurement for those hours. The vertical axis is drawn to fit the middle 99%
						of the hours, so the most extreme few are off the scale. Drag the ends of the bar under the figure
						to zoom in.
					</p>
				{:else if view && averaged}
					<p class="muted caption">
						Each {isolating ? 'dark grey' : 'grey'} point is one {AVERAGE_NOUN[view.settings.averaging]} average{isolating
							? ' of the isolated signal'
							: ''}.
						{#if measuring && run}
							The <span class="key key--trend">red line</span> is the trend AQEval fitted.
							<span class="key key--break">Blue</span> marks where it changes: a dashed line where one long
							stretch gives way to the next, and a shaded band where the change was over within weeks. The
							black lines mark the events. The break points from step 4 are switched off: click
							<em>Break points</em> in the key to see them.
						{:else if run}
							The <span class="key key--breakpoint">dashed lines</span> are the break points, each inside a
							shaded band: the range of dates it could fall in. The solid black lines mark the events.
						{/if}
						The vertical axis is drawn to fit the middle 99% of the points, so the most extreme few are off the
						scale. Drag the bar under the figure to zoom{run ? ', and click a name in the key to show or hide it' : ''}.
					</p>
					{#if isolating}
						<label class="check">
							<input type="checkbox" bind:checked={showMeasured} />
							Show the measurements underneath, in lighter grey
						</label>
						{#if view.runs?.formula}
							<p class="muted model">
								Model removed from the measurements: <code>{view.runs.formula}</code>
							</p>
						{/if}
					{/if}
				{/if}
				{#if failure}
					<p class="alert alert--error" role="alert">
						The results for this set-up could not be loaded{view ? ', so the figure still shows the last one' : ''}.
						Check your connection, or try another setting. ({failure})
					</p>
				{/if}

				{#if experiment.nextStage}
					{@const next = experiment.nextStage}
					<div class="next">
						<button class="btn" onclick={() => experiment.goTo(next)}>
							Next: {STAGE_ACTION[next].toLowerCase()} <span aria-hidden="true">→</span>
						</button>
						<p class="muted">{STAGE_TEASER[next]}</p>
					</div>
				{/if}
			</section>

			{#if !wide}
				<section class="card" aria-label="Code">
					<CodePanel />
				</section>
			{/if}

			<section class="card" aria-labelledby="found-title">
				<h2 id="found-title" class="card__title">{NOTES_TITLE[view?.stage ?? experiment.stage]}</h2>
				{#if view && measuring && run}
					<Report
						{run}
						events={manifest.events}
						{pollutant}
						maxBreaks={manifest.maxBreaksToQuantify}
						language={view.language}
					/>
				{:else if view && !measuring}
					<StageNotes
						stage={view.stage}
						settings={view.settings}
						hourly={view.hourly}
						{averaged}
						measured={measuredAveraged}
						{isolating}
						{run}
						events={manifest.events}
						{pollutant}
						{siteName}
						years={manifest.years}
						language={view.language}
					/>
				{:else if failure}
					<p class="muted">Nothing to report until a result loads.</p>
				{:else}
					<div class="skeleton" style:height="12rem"></div>
				{/if}
			</section>

			{#if measuring || pins.length}
				<section class="card" aria-labelledby="compare-title">
					<h2 id="compare-title" class="card__title">Compare results</h2>
					<Compare
						{pins}
						{siteName}
						onLoad={(p) => {
							experiment.language = p.language;
							experiment.set(p.settings);
							experiment.goTo('changes');
						}}
						onRemove={(p) => (pins = pins.filter((other) => other.id !== p.id))}
					/>
				</section>
			{/if}

			<section class="card" aria-labelledby="sites-title">
				<div class="card__head">
					<h2 id="sites-title" class="card__title">Compare the sites</h2>
					<Segmented
						legend="Label each site with"
						hideLegend
						small
						name="metric"
						options={[
							{ value: 'before', label: 'Year before' },
							{ value: 'after', label: 'Year after' },
							{ value: 'change', label: 'Change' }
						]}
						value={metric}
						onChange={(value) => (metric = value)}
					/>
				</div>
				<p class="muted">
					{#if metric === 'change'}
						Change in average {POLLUTANT_LABEL[experiment.settings.pollutant]} from the twelve months before
						the Clean Air Zone started{caz ? ` (${formatDate(caz.date)})` : ''} to the twelve months after.
					{:else}
						Average {POLLUTANT_LABEL[experiment.settings.pollutant]} over the twelve months
						{metric} the Clean Air Zone started{caz ? ` (${formatDate(caz.date)})` : ''}.
					{/if}
					These are plain averages of the measurements, with nothing removed. Click a site to analyse it.
				</p>
				<SiteMap
					sites={manifest.sites}
					{boundary}
					pollutant={experiment.settings.pollutant}
					selected={experiment.settings.site}
					background={experiment.settings.background}
					{metric}
					onSelect={(code) =>
						experiment.set({
							site: code,
							...(experiment.settings.background === code ? { background: null } : {})
						})}
				/>

				<h3 class="sub">Average for each year</h3>
				{#if annual}
					<Chart
						option={annual}
						height="230px"
						description={`Annual mean ${POLLUTANT_LABEL[experiment.settings.pollutant]} at each site, ${manifest.years[0]} to ${manifest.years[1]}`}
					/>
				{/if}
				{#if selectedSite}
					<p class="muted capture">
						Hours with a valid {POLLUTANT_LABEL[experiment.settings.pollutant]} measurement at
						{selectedSite.name}:
						{#each captureYears as [year, share], index (year)}
							<span class:low={(share ?? 0) < 75}>{year} {share ?? 0}%</span>{index < captureYears.length - 1
								? ', '
								: '.'}
						{/each}
						A year under 75% is usually treated as incomplete.
					</p>
				{/if}

				<h3 class="sub">NO₂ as a share of NOx, each month</h3>
				{#if ratio}
					<Chart
						option={ratio}
						height="230px"
						description={`Monthly ratio of NO₂ to NOx at each site, ${manifest.years[0]} to ${manifest.years[1]}`}
					/>
				{/if}
				<p class="muted capture">
					The measurements only, with no analysis: each month's average NO₂ divided by its average NOx, over
					the hours when both were measured. A month under 75% complete is left out. NOx is NO plus NO₂.
					Next to a busy road most of it is NO, fresh from exhausts, so NO₂'s share is low. Away from the
					traffic the NO has had time to react with ozone and become NO₂, so the share is higher.
				</p>
			</section>
		</div>

		{#if wide}
			<aside class="bench__code card" aria-label="Code">
				<CodePanel />
			</aside>
		{/if}
	</div>

	<section class="band" id="tasks" aria-labelledby="tasks-title">
		<h2 id="tasks-title">Tasks</h2>
		<p class="band__lede">
			Each task can be answered on this page by moving the dials. Most need every step of the analysis
			added. Then try it in code.
		</p>
		<Tasks />
	</section>

	<section class="band about" id="about" aria-labelledby="about-title">
		<h2 id="about-title">About this page</h2>
		<div class="about__columns">
			<div>
				<h3>Where the results come from</h3>
				<p>
					Nothing is calculated in your browser. Every set-up the dials allow was run ahead of time, and
					the page looks up the answer. The code panel shows what was run.
				</p>
				<p>
					There are two sets of answers, one computed in R and one in Python, and the page shows the set
					for the language chosen in the code panel. The two languages do not give identical answers,
					for one reason: they fit the model that isolates the signal with different libraries. The two
					fits are close, about 1 µg/m³ apart, and that is enough to move where a change is found. Switch
					language with a result pinned to see how much it matters. With nothing isolated the two agree.
				</p>
				<p class="muted">
					R: {manifest.engines.r}. Python: {manifest.engines.python}. One step of the R set was not run in
					R: the search for break points, which takes R up to a quarter of an hour a time. It was run with
					the Python version, which returns the same break points.
				</p>
			</div>
			<div>
				<h3>Run it yourself</h3>
				<p>
					<a href={CODESPACE_URL} target="_blank" rel="noopener">Open the workspace</a> to get RStudio and a
					Python notebook in your browser, with the packages already installed. It needs a free GitHub
					account. The R script and the notebook there are the same steps you see on this page.
				</p>
				<p>
					The original step-by-step walkthrough is in the
					<a href={REPOSITORY_URL} target="_blank" rel="noopener">repository</a>.
				</p>
			</div>
			<div>
				<h3>Data and credits</h3>
				<p>
					Measurements: Automatic Urban and Rural Network, © Crown copyright Defra, via
					<a href="https://uk-air.defra.gov.uk" target="_blank" rel="noopener">uk-air.defra.gov.uk</a>,
					Open Government Licence. Clean Air Zone boundary: City of Bradford Metropolitan District Council,
					Open Government Licence v3.0. Results generated {formatDate(manifest.generated)}.
				</p>
				<p>
					Method: Ropkins, Walker and Tate, <em>AQEval</em>.
					<a href="https://doi.org/10.21105/joss.08839" target="_blank" rel="noopener"
						>Journal of Open Source Software</a
					>.
				</p>
			</div>
		</div>
	</section>
{/if}

<style>
	.intro {
		max-width: 62rem;
		padding: var(--space-5) var(--space-4) var(--space-2);
	}

	.intro h1 {
		font-size: clamp(1.5rem, 1.1rem + 1.6vw, 2.25rem);
		margin-bottom: var(--space-3);
	}

	.intro__lede {
		font-size: 1.0625rem;
	}

	.intro__how {
		margin-bottom: 0;
	}

	.load-error {
		margin: var(--space-4);
	}

	/* One column on a phone; dials beside results on a laptop; and on a wide
	   screen the code sits alongside too, so a dial, the figure and the code
	   are all in view at once. */
	.bench {
		display: grid;
		grid-template-columns: minmax(0, 1fr);
		gap: var(--space-4);
		padding: var(--space-4);
		align-items: start;
	}

	.bench__main {
		display: grid;
		/* Without an explicit column a grid grows to fit its widest child, and
		   a chart or a table would push the page sideways. */
		grid-template-columns: minmax(0, 1fr);
		gap: var(--space-4);
		min-width: 0;
	}

	.bench__code {
		min-width: 0;
	}

	@media (min-width: 900px) {
		.bench {
			grid-template-columns: 300px minmax(0, 1fr);
		}

		.bench__controls {
			position: sticky;
			top: calc(var(--header-height) + var(--space-4));
			max-height: calc(100vh - var(--header-height) - var(--space-6));
			overflow-y: auto;
		}

	}

	@media (min-width: 1360px) {
		.bench {
			grid-template-columns: 300px minmax(0, 1fr) 480px;
		}

		.bench__code {
			position: sticky;
			top: calc(var(--header-height) + var(--space-4));
			max-height: calc(100vh - var(--header-height) - var(--space-6));
			overflow-y: auto;
		}
	}

	.card__head {
		display: flex;
		flex-wrap: wrap;
		align-items: flex-start;
		justify-content: space-between;
		gap: var(--space-2) var(--space-3);
		margin-bottom: var(--space-2);
	}

	.card__head .card__title {
		margin-bottom: 0;
	}

	/* The title gives way before the control beside it is pushed underneath. */
	.card__head > :first-child {
		flex: 1 1 12rem;
		min-width: 0;
	}

	.setup {
		margin: 2px 0 0;
	}

	.figure {
		transition: opacity 0.15s ease;
	}

	.figure--loading {
		opacity: 0.5;
	}

	.caption {
		margin: var(--space-2) 0;
	}

	.key {
		font-weight: 700;
	}

	.key--trend {
		color: var(--colour-trend);
	}

	.key--break {
		color: var(--colour-break);
	}

	.key--breakpoint {
		color: #505a5f;
	}

	/* The way on to the next step, at the foot of the figure. */
	.next {
		display: flex;
		flex-wrap: wrap;
		align-items: center;
		gap: var(--space-2) var(--space-3);
		margin-top: var(--space-3);
		padding-top: var(--space-3);
		border-top: 1px solid var(--border-subtle);
	}

	.next p {
		flex: 1 1 16rem;
		margin: 0;
		font-size: 0.875rem;
	}

	.check {
		display: flex;
		align-items: center;
		gap: var(--space-2);
		font-size: 0.875rem;
		cursor: pointer;
		width: fit-content;
	}

	.check input {
		width: 1.125rem;
		height: 1.125rem;
		accent-color: var(--text);
	}

	.model {
		margin: var(--space-2) 0 0;
	}

	.sub {
		font-size: 0.9375rem;
		margin: var(--space-4) 0 var(--space-1);
	}

	.capture {
		margin: var(--space-2) 0 0;
	}

	.low {
		color: var(--colour-rise);
		font-weight: 700;
	}

	.band {
		padding: var(--space-5) var(--space-4);
	}

	.band > h2 {
		font-size: 1.5rem;
		margin-bottom: var(--space-2);
	}

	.band__lede {
		max-width: 50rem;
	}

	.about {
		background: var(--surface);
		border-top: 1px solid var(--border-subtle);
	}

	.about__columns {
		display: grid;
		grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
		gap: var(--space-5);
		font-size: 0.9375rem;
	}

	.about h3 {
		font-size: 1rem;
		margin-bottom: var(--space-2);
	}
</style>
