<script lang="ts" generics="T extends string | number">
	/**
	 * A row of mutually exclusive choices, for when there are few enough to
	 * show them all. Radio buttons underneath, so the keyboard and screen
	 * readers get the behaviour they expect.
	 */
	interface Option {
		value: T;
		label: string;
	}

	interface Props {
		legend: string;
		/** Hide the legend visually when a heading already names the group. */
		hideLegend?: boolean;
		name: string;
		options: Option[];
		value: T;
		onChange: (value: T) => void;
		small?: boolean;
	}

	let { legend, hideLegend = false, name, options, value, onChange, small = false }: Props = $props();
</script>

<fieldset class="segmented" class:segmented--small={small}>
	<legend class:visually-hidden={hideLegend}>{legend}</legend>
	<div class="segmented__row">
		{#each options as option (option.value)}
			<label class="segmented__option" class:segmented__option--on={option.value === value}>
				<input
					type="radio"
					{name}
					class="visually-hidden"
					checked={option.value === value}
					onchange={() => onChange(option.value)}
				/>
				{option.label}
			</label>
		{/each}
	</div>
</fieldset>

<style>
	.segmented {
		border: 0;
		margin: 0;
		padding: 0;
		min-width: 0;
	}

	legend {
		padding: 0;
		margin-bottom: var(--space-1);
		font-size: 0.875rem;
		font-weight: 600;
	}

	.segmented__row {
		display: flex;
		border: 2px solid var(--text);
		border-radius: var(--radius);
		overflow: hidden;
		background: var(--surface);
	}

	.segmented__option {
		flex: 1;
		padding: var(--space-2) var(--space-2);
		text-align: center;
		font-weight: 600;
		font-size: 0.9375rem;
		cursor: pointer;
		white-space: nowrap;
		border-left: 1px solid var(--border);
	}

	.segmented__option:first-child {
		border-left: 0;
	}

	.segmented__option:hover {
		background: var(--surface-sunken);
	}

	.segmented__option--on,
	.segmented__option--on:hover {
		background: var(--text);
		color: var(--text-inverse);
	}

	/* The input is hidden, so the label has to show the focus ring for it. */
	.segmented__option:has(:focus-visible) {
		outline: 3px solid var(--focus);
		outline-offset: -3px;
	}

	.segmented--small .segmented__row {
		border-width: 1px;
	}

	.segmented--small .segmented__option {
		padding: 2px var(--space-2);
		font-size: 0.8125rem;
	}
</style>
