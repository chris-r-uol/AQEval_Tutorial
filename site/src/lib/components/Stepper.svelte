<script lang="ts">
	/**
	 * Where the reader is in the analysis.
	 *
	 * One button per step, in order. The steps added so far are filled in;
	 * choosing an earlier one takes the later ones away again, and choosing a
	 * later one adds everything up to it.
	 */
	import { STAGES } from '$lib/code';
	import { experiment } from '$lib/experiment.svelte';
	import { STAGE_LABEL } from '$lib/stages';

	let current = $derived(STAGES.indexOf(experiment.stage));
</script>

<nav class="stepper" aria-label="Steps of the analysis">
	<ol>
		{#each STAGES as stage, index (stage)}
			<li class:done={index < current} class:current={index === current}>
				<button
					type="button"
					aria-current={index === current ? 'step' : undefined}
					onclick={() => experiment.goTo(stage)}
				>
					<span class="number" aria-hidden="true">{index + 1}</span>
					<span class="label">{STAGE_LABEL[stage]}</span>
				</button>
			</li>
		{/each}
	</ol>
</nav>

<style>
	.stepper {
		margin: 0 0 var(--space-3);
		overflow-x: auto;
	}

	ol {
		display: flex;
		margin: 0;
		padding: 0;
		list-style: none;
		min-width: max-content;
	}

	li {
		flex: 1 1 0;
		position: relative;
	}

	/* The line joining each step to the one before. */
	li + li::before {
		content: '';
		position: absolute;
		top: calc(var(--space-1) + 0.75rem);
		right: 50%;
		width: 100%;
		height: 2px;
		background: var(--border-subtle);
		z-index: 0;
	}

	li.done + li::before,
	li.done + li.current::before {
		background: var(--text);
	}

	button {
		position: relative;
		z-index: 1;
		display: flex;
		flex-direction: column;
		align-items: center;
		gap: 2px;
		width: 100%;
		padding: var(--space-1) var(--space-2);
		border: 0;
		background: none;
		color: var(--text-muted);
		font: inherit;
		font-size: 0.8125rem;
		cursor: pointer;
	}

	button:hover .label {
		text-decoration: underline;
	}

	button:focus-visible {
		outline: 3px solid var(--focus);
		outline-offset: 0;
	}

	.number {
		display: grid;
		place-items: center;
		width: 1.5rem;
		height: 1.5rem;
		border-radius: 50%;
		border: 2px solid var(--border);
		background: var(--surface);
		font-weight: 700;
		font-size: 0.75rem;
	}

	.done .number,
	.current .number {
		border-color: var(--text);
		background: var(--text);
		color: var(--text-inverse);
	}

	.done button,
	.current button {
		color: var(--text);
	}

	.current .label {
		font-weight: 700;
	}

	.current .number {
		box-shadow: 0 0 0 3px var(--dial-strong);
	}

	/* On a phone the five names do not fit side by side: name only the step
	   the reader is on. The others keep their names for screen readers. */
	@media (max-width: 600px) {
		ol {
			min-width: 0;
		}

		li:not(.current) .label {
			position: absolute;
			width: 1px;
			height: 1px;
			overflow: hidden;
			clip-path: inset(50%);
			white-space: nowrap;
		}

		li.current {
			flex: 3 1 auto;
		}

		.current .label {
			white-space: nowrap;
		}
	}
</style>
