<script lang="ts">
	/**
	 * The code behind the figure.
	 *
	 * Shows the R or Python that produces exactly what is on screen. Every
	 * value a dial controls is highlighted, and flashes when its dial moves, so
	 * turning a dial teaches which argument it is.
	 */
	import { untrack } from 'svelte';
	import { plain, script, segments, steps, type Settings } from '$lib/code';
	import { CODESPACE_URL } from '$lib/config';
	import { experiment } from '$lib/experiment.svelte';
	import { highlight } from '$lib/highlight';
	import Segmented from './Segmented.svelte';

	let settings = $derived(experiment.settings);
	let language = $derived(experiment.language);
	let list = $derived(steps(settings, language));

	// Below h = 0.3 an 8-hour series is slow in both languages: R's search for
	// break points slows sharply, and in either language measuring k breaks
	// takes 9^k model fits.
	let slow = $derived(settings.averaging === '8 hour' && settings.h < 0.3);

	// Which dials changed most recently, so their values can flash.
	let changed = $state<string[]>([]);
	let previous: Settings | null = null;
	let timer: ReturnType<typeof setTimeout> | undefined;

	$effect(() => {
		const current = settings;
		untrack(() => {
			if (previous) {
				const keys = (Object.keys(current) as Array<keyof Settings>).filter(
					(key) => current[key] !== previous![key]
				);
				if (keys.length) {
					changed = keys;
					clearTimeout(timer);
					timer = setTimeout(() => (changed = []), 1400);
				}
			}
			previous = current;
		});
	});

	/** One step's code as HTML: syntax colours, with dial values marked. */
	function render(code: string, flashing: string[]): string {
		return segments(code)
			.map((part) =>
				part.dial
					? `<mark class="dial${flashing.includes(part.dial) ? ' dial--flash' : ''}">${highlight(part.text)}</mark>`
					: highlight(part.text)
			)
			.join('');
	}

	let copied = $state(false);

	async function copy() {
		await navigator.clipboard.writeText(plain(script(settings, language)));
		copied = true;
		setTimeout(() => (copied = false), 2000);
	}
</script>

<div class="code">
	<div class="code__head">
		<h2 class="code__title">The code for this figure</h2>
		<Segmented
			legend="Language"
			hideLegend
			small
			name="language"
			options={[
				{ value: 'r', label: 'R' },
				{ value: 'python', label: 'Python' }
			]}
			value={language}
			onChange={(value) => (experiment.language = value)}
		/>
	</div>

	<p class="muted code__hint">
		<mark class="dial">Highlighted</mark> values come from the dials. Change a dial and watch the code change.
	</p>

	<!-- The steps are rendered from our own templates, with values taken only
	     from the options the dials offer, so the markup is safe to insert. -->
	<!-- Focusable so the keyboard can scroll it when a line is wider than the panel. -->
	<!-- svelte-ignore a11y_no_noninteractive_tabindex -->
	<pre class="code__block" tabindex="0" aria-label={`${language === 'r' ? 'R' : 'Python'} code`}><code
			>{#each list as step, index (step.id)}<span class="tok-comment"># {index + 1}. {step.title}</span>
{@html render(step.code, changed)}{index < list.length - 1 ? '\n\n' : ''}{/each}</code
		></pre>

	<div class="code__actions">
		<button class="btn btn--secondary btn--small" onclick={copy}>{copied ? 'Copied' : 'Copy the code'}</button>
		<a class="btn btn--small" href={CODESPACE_URL} target="_blank" rel="noopener">Run it yourself</a>
	</div>
	<p class="muted code__note">
		<strong>Run it yourself</strong> opens a ready-made workspace in your browser with
		{language === 'r' ? 'RStudio' : 'a Python notebook'} and every package installed. Nothing to set up.
	</p>
	{#if slow}
		<p class="muted code__note">
			<strong>This set-up is slow to run yourself.</strong> With 8-hour data and <code>h</code> below 0.3,
			{language === 'r'
				? "R's search for break points takes tens of minutes"
				: 'measuring the changes can take several minutes'}. Daily averages take a minute or two.
		</p>
	{/if}
	<span class="visually-hidden" role="status">{copied ? 'Code copied' : ''}</span>
</div>

<style>
	.code__head {
		display: flex;
		align-items: center;
		justify-content: space-between;
		gap: var(--space-3);
		margin-bottom: var(--space-2);
	}

	.code__title {
		font-size: 1rem;
	}

	.code__hint {
		margin-bottom: var(--space-2);
	}

	.code__block {
		margin: 0;
		padding: var(--space-3);
		background: var(--surface-code);
		border: 1px solid var(--border-subtle);
		border-radius: var(--radius);
		font-family: var(--font-mono);
		font-size: 0.75rem;
		line-height: 1.65;
		overflow-x: auto;
		white-space: pre;
		tab-size: 2;
	}

	.code__block code {
		font: inherit;
		background: none;
		padding: 0;
	}

	.code__actions {
		display: flex;
		flex-wrap: wrap;
		gap: var(--space-2);
		margin: var(--space-3) 0 var(--space-2);
	}

	.code__note {
		margin: 0 0 var(--space-2);
	}

	mark.dial,
	.code__block :global(mark.dial) {
		background: var(--dial);
		color: inherit;
		border-radius: 3px;
		padding: 1px 2px;
		box-shadow: inset 0 -2px 0 var(--dial-strong);
	}

	.code__block :global(mark.dial--flash) {
		animation: flash 1.4s ease-out;
	}

	@keyframes flash {
		0%,
		35% {
			background: var(--dial-strong);
			box-shadow: 0 0 0 3px var(--dial-strong);
		}
		100% {
			background: var(--dial);
			box-shadow: inset 0 -2px 0 var(--dial-strong);
		}
	}

	.code__block :global(.tok-comment) {
		color: #505a5f;
		font-style: italic;
	}

	.code__block :global(.tok-string) {
		color: #00703c;
	}

	.code__block :global(.tok-number) {
		color: #b0420a;
	}

	.code__block :global(.tok-keyword) {
		color: #4c2c92;
		font-weight: 600;
	}

	.code__block :global(.tok-call) {
		color: #003078;
	}
</style>
