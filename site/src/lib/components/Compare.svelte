<script lang="ts">
	/**
	 * Pinned results, side by side.
	 *
	 * An experiment is a comparison: the same site with and without the weather
	 * removed, a roadside against a background site, one sensitivity against
	 * another. Pinning keeps a result on the figure and in this table while the
	 * dials move on.
	 */
	import { summarise } from '$lib/analysis';
	import { AVERAGING_LABEL, POLLUTANT_LABEL, describeIsolation, formatNumber, formatSigned } from '$lib/format';
	import type { Pin } from '$lib/types';

	interface Props {
		pins: Pin[];
		siteName: (code: string) => string;
		onLoad: (pin: Pin) => void;
		onRemove: (pin: Pin) => void;
	}

	let { pins, siteName, onLoad, onRemove }: Props = $props();
</script>

{#if pins.length === 0}
	<p class="muted empty">
		Nothing pinned yet. Use <strong>Pin this result</strong> above the figure, then change a dial: the pinned
		trend stays on the figure as a dashed line, and both results are listed here.
	</p>
{:else}
	<div class="table-scroll">
		<table>
			<thead>
				<tr>
					<th scope="col">Pin</th>
					<th scope="col">Site</th>
					<th scope="col">Set-up</th>
					<th scope="col" class="num">Changes</th>
					<th scope="col" class="num">Largest change</th>
					<th scope="col" class="num">Start to end</th>
					<th scope="col"><span class="visually-hidden">Actions</span></th>
				</tr>
			</thead>
			<tbody>
				{#each pins as pin (pin.id)}
					{@const summary = summarise(pin.run)}
					<tr>
						<th scope="row">
							<span class="swatch" style:border-top-color={pin.colour} aria-hidden="true"></span>
							{pin.label}
						</th>
						<td>
							{siteName(pin.settings.site)}
							<span class="muted">{POLLUTANT_LABEL[pin.settings.pollutant]}</span>
						</td>
						<td class="setup">
							{describeIsolation(pin.settings)}
							<span class="muted">
								· {AVERAGING_LABEL[pin.settings.averaging]} · h {pin.settings.h} ·
								{pin.language === 'r' ? 'R' : 'Python'}
							</span>
						</td>
						<td class="num">{summary.changes}</td>
						<td class="num">
							{summary.largest?.material ? formatSigned(summary.largest.percent, 0, '%') : '–'}
						</td>
						<td class="num">
							{formatSigned(summary.overallPercent, 0, '%')}
							<span class="muted">({formatNumber(summary.start)} → {formatNumber(summary.end)})</span>
						</td>
						<td class="actions">
							<button class="btn btn--secondary btn--small" onclick={() => onLoad(pin)}>
								Show<span class="visually-hidden"> pin {pin.label} on the dials</span>
							</button>
							<button class="btn btn--secondary btn--small" onclick={() => onRemove(pin)}>
								Remove<span class="visually-hidden"> pin {pin.label}</span>
							</button>
						</td>
					</tr>
				{/each}
			</tbody>
		</table>
	</div>
{/if}

<style>
	.empty {
		margin: 0;
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
		vertical-align: middle;
	}

	thead th {
		border-bottom: 2px solid var(--text);
		white-space: nowrap;
	}

	.num {
		text-align: right;
		white-space: nowrap;
	}

	.setup {
		min-width: 12rem;
	}

	.swatch {
		display: inline-block;
		width: 1.25rem;
		border-top: 3px dashed;
		vertical-align: middle;
		margin-right: var(--space-1);
	}

	.actions {
		white-space: nowrap;
		text-align: right;
	}

	.actions .btn + .btn {
		margin-left: var(--space-1);
	}
</style>
