<script lang="ts">
	/**
	 * The dials. Each group is one step of the analysis, in the order the code
	 * runs them, so setting up the experiment reads top to bottom like the
	 * script it produces.
	 */
	import { experiment } from '$lib/experiment.svelte';
	import { AVERAGING_LABEL, POLLUTANT_LABEL, POLLUTANT_NAME, formatDuration } from '$lib/format';
	import { SITE_COLOURS } from '$lib/charts';
	import type { Manifest } from '$lib/types';
	import Help from './Help.svelte';
	import Segmented from './Segmented.svelte';

	interface Props {
		manifest: Manifest;
	}

	let { manifest }: Props = $props();

	let settings = $derived(experiment.settings);
	let pollutant = $derived(POLLUTANT_LABEL[settings.pollutant]);

	let siteName = (code: string) => manifest.sites.find((site) => site.code === code)?.name ?? code;

	// The background is one <select>, but its values are null, a column name or
	// a site code; "none" stands in for null in the markup.
	let backgroundOptions = $derived([
		{ value: 'none', label: 'Nothing' },
		{ value: 'air_temp', label: 'Air temperature' },
		...manifest.controlSites
			.filter((code) => code !== settings.site)
			.map((code) => ({ value: code, label: `${pollutant} at ${siteName(code)}` }))
	]);

	function chooseSite(code: string) {
		// A site cannot be its own background.
		experiment.set({ site: code, ...(settings.background === code ? { background: null } : {}) });
	}

	// The slider moves through the precomputed values of h by position.
	let hIndex = $derived(manifest.h.indexOf(settings.h));
	let spanDays = $derived((manifest.years[1] - manifest.years[0] + 1) * 365.25);
</script>

<div class="controls">
	<h2 class="controls__title">Set up the experiment</h2>

	<ol class="steps">
		<li class="step">
			<h3><span class="step__number">1</span> Choose a site</h3>
			<div class="sites" role="radiogroup" aria-label="Monitoring site">
				{#each manifest.sites as site (site.code)}
					<label class="site" class:site--on={site.code === settings.site}>
						<input
							type="radio"
							name="site"
							class="visually-hidden"
							checked={site.code === settings.site}
							onchange={() => chooseSite(site.code)}
						/>
						<span class="site__dot" style:background={SITE_COLOURS[site.type]} aria-hidden="true"></span>
						<span class="site__text">
							<span class="site__name">{site.name}</span>
							<span class="site__meta">
								{site.type === 'Urban Traffic' ? 'Roadside' : 'Background'} ·
								{site.insideCaz ? 'inside the Clean Air Zone' : `${site.kmFromBradford.toFixed(0)} km away, no zone`}
							</span>
						</span>
						<code class="site__code">{site.code.toLowerCase()}</code>
					</label>
				{/each}
			</div>
			<Help>
				<p>
					Each site is a station on the national monitoring network (AURN). A <strong>roadside</strong> site
					sits next to traffic, so it responds to what vehicles are doing. A <strong>background</strong> site
					sits away from busy roads and shows the air the whole area shares.
				</p>
				<p>
					Bradford has one station with a long record, on Mayo Avenue. The others are its nearest
					neighbours, none of them in a charging zone, so they show what happened where there was no
					Clean Air Zone.
				</p>
			</Help>
		</li>

		<li class="step">
			<h3><span class="step__number">2</span> Choose a pollutant</h3>
			<Segmented
				legend="Pollutant"
				hideLegend
				name="pollutant"
				options={manifest.pollutants.map((value) => ({ value, label: POLLUTANT_LABEL[value] }))}
				value={settings.pollutant}
				onChange={(value) => experiment.set({ pollutant: value })}
			/>
			<Help>
				<p>
					<strong>NO₂</strong> ({POLLUTANT_NAME.no2}) is the pollutant with a legal limit, and the one Clean
					Air Zones exist to bring down. <strong>NOx</strong> ({POLLUTANT_NAME.nox}) is NO₂ plus nitric oxide,
					and follows exhaust emissions more directly.
				</p>
			</Help>
		</li>

		<li class="step">
			<h3><span class="step__number">3</span> Isolate the signal</h3>
			<label class="check">
				<input
					type="checkbox"
					checked={settings.deseason}
					onchange={(event) => experiment.set({ deseason: event.currentTarget.checked })}
				/>
				Remove the daily and yearly cycle
			</label>
			<label class="check">
				<input
					type="checkbox"
					checked={settings.deweather}
					onchange={(event) => experiment.set({ deweather: event.currentTarget.checked })}
				/>
				Remove the effect of the wind
			</label>
			<label class="select">
				<span>Also account for</span>
				<select
					value={settings.background ?? 'none'}
					onchange={(event) => {
						const value = event.currentTarget.value;
						experiment.set({ background: value === 'none' ? null : value });
					}}
				>
					{#each backgroundOptions as option (option.value)}
						<option value={option.value}>{option.label}</option>
					{/each}
				</select>
			</label>
			<Help>
				<p>
					Pollution is higher at rush hour, in winter and on still days. None of that is a change in
					emissions. AQEval fits a model of those patterns and subtracts it, so what is left is the part
					a policy could have caused.
				</p>
				<p>
					<strong>Also account for</strong> adds one more thing to the model. Choosing a background site
					removes whatever the two sites have in common, such as region-wide changes, and leaves what is
					particular to this site.
				</p>
				<p>Untick everything to analyse the measurements as they are.</p>
			</Help>
		</li>

		<li class="step">
			<h3><span class="step__number">4</span> Average the data</h3>
			<Segmented
				legend="Averaging period"
				hideLegend
				name="averaging"
				options={manifest.averaging.map((value) => ({ value, label: AVERAGING_LABEL[value] }))}
				value={settings.averaging}
				onChange={(value) => experiment.set({ averaging: value })}
			/>
			<Help>
				<p>
					The data arrives as one value an hour. Averaging into longer blocks smooths out short-lived
					spikes. Shorter blocks keep more detail and can date a change more exactly, but are noisier and
					slower to analyse.
				</p>
			</Help>
		</li>

		<li class="step">
			<h3><span class="step__number">5</span> Set the sensitivity</h3>
			<label class="slider">
				<span class="slider__label">
					<code>h = {settings.h}</code>
					<span class="muted">breaks at least {formatDuration(settings.h * spanDays)} apart</span>
				</span>
				<input
					type="range"
					min="0"
					max={manifest.h.length - 1}
					step="1"
					value={hIndex}
					list="h-values"
					aria-valuetext={`h = ${settings.h}`}
					oninput={(event) => experiment.set({ h: manifest.h[Number(event.currentTarget.value)] })}
				/>
				<span class="slider__ends" aria-hidden="true">
					<span>More sensitive</span>
					<span>Less sensitive</span>
				</span>
			</label>
			<datalist id="h-values">
				{#each manifest.h as _, index (index)}
					<option value={index}></option>
				{/each}
			</datalist>
			<Help>
				<p>
					<code>h</code> is the shortest stretch allowed between two breaks, as a fraction of the whole
					series. With {manifest.years[1] - manifest.years[0] + 1} years of data, <code>h = 0.3</code>
					means no two breaks closer than about {formatDuration(0.3 * spanDays)}.
				</p>
				<p>
					A smaller <code>h</code> is more sensitive: it can find more changes, closer together. It can also
					split one real change into several, and it takes much longer to run.
				</p>
				<p>
					The dial starts at 0.3, a quick first look. For 8-hour data we would normally use
					<code>h = 0.12</code>, the most sensitive setting here. Work down towards it and watch what
					appears.
				</p>
			</Help>
		</li>
	</ol>

	<button class="btn btn--secondary" onclick={() => experiment.reset()} disabled={experiment.isDefault}>
		Reset to the starting set-up
	</button>
</div>

<style>
	.controls__title {
		font-size: 1.125rem;
		margin-bottom: var(--space-3);
	}

	.steps {
		list-style: none;
		margin: 0 0 var(--space-4);
		padding: 0;
	}

	.step {
		padding: var(--space-3) 0;
		border-top: 1px solid var(--border-subtle);
	}

	.step h3 {
		font-size: 0.9375rem;
		margin-bottom: var(--space-2);
		display: flex;
		align-items: center;
		gap: var(--space-2);
	}

	.step__number {
		display: inline-grid;
		place-items: center;
		width: 1.5rem;
		height: 1.5rem;
		border-radius: 50%;
		background: var(--text);
		color: var(--text-inverse);
		font-size: 0.8125rem;
	}

	.sites {
		display: grid;
		gap: var(--space-1);
	}

	.site {
		display: flex;
		align-items: center;
		gap: var(--space-2);
		padding: var(--space-2);
		border: 2px solid var(--border-subtle);
		border-radius: var(--radius);
		cursor: pointer;
	}

	.site:hover {
		border-color: var(--border);
	}

	.site--on,
	.site--on:hover {
		border-color: var(--text);
		background: var(--surface-sunken);
	}

	.site:has(:focus-visible) {
		outline: 3px solid var(--focus);
	}

	.site__dot {
		width: 0.75rem;
		height: 0.75rem;
		border-radius: 50%;
		flex-shrink: 0;
	}

	.site__text {
		flex: 1;
		min-width: 0;
		line-height: 1.25;
	}

	.site__name {
		display: block;
		font-weight: 600;
		font-size: 0.9375rem;
	}

	.site__meta {
		display: block;
		color: var(--text-muted);
		font-size: 0.8125rem;
	}

	.site__code {
		flex-shrink: 0;
	}

	.check {
		display: flex;
		align-items: center;
		gap: var(--space-2);
		padding: var(--space-1) 0;
		font-size: 0.9375rem;
		cursor: pointer;
	}

	.check input {
		width: 1.25rem;
		height: 1.25rem;
		accent-color: var(--text);
		flex-shrink: 0;
	}

	.select {
		display: block;
		margin-top: var(--space-2);
		font-size: 0.9375rem;
	}

	.select span {
		display: block;
		font-weight: 600;
		font-size: 0.875rem;
		margin-bottom: var(--space-1);
	}

	.select select {
		width: 100%;
		padding: var(--space-2);
		border: 2px solid var(--text);
		border-radius: 0;
		background: var(--surface);
		color: var(--text);
	}

	.slider {
		display: block;
	}

	.slider__label {
		display: flex;
		align-items: baseline;
		justify-content: space-between;
		gap: var(--space-2);
	}

	.slider__label code {
		font-size: 1rem;
		font-weight: 700;
	}

	.slider input {
		width: 100%;
		margin: var(--space-2) 0 0;
		accent-color: var(--text);
	}

	.slider__ends {
		display: flex;
		justify-content: space-between;
		color: var(--text-muted);
		font-size: 0.75rem;
	}
</style>
