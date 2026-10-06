<script lang="ts">
	/**
	 * Thin ECharts wrapper.
	 *
	 * Only the pieces of ECharts actually used are imported, so the bundle
	 * carries a few hundred KB rather than the full library.
	 */
	import { onMount } from 'svelte';
	import * as echarts from 'echarts/core';
	import { LineChart, ScatterChart } from 'echarts/charts';
	import {
		GridComponent,
		TooltipComponent,
		LegendComponent,
		MarkAreaComponent,
		MarkLineComponent,
		DataZoomComponent
	} from 'echarts/components';
	import { CanvasRenderer } from 'echarts/renderers';
	import type { EChartsOption } from 'echarts';

	// Everything a chart on this page uses must be registered here. ECharts
	// throws on an unregistered series type or component rather than ignoring it.
	echarts.use([
		LineChart,
		ScatterChart,
		GridComponent,
		TooltipComponent,
		LegendComponent,
		MarkAreaComponent,
		MarkLineComponent,
		DataZoomComponent,
		CanvasRenderer
	]);

	interface Props {
		option: EChartsOption;
		height?: string;
		/** Describes the chart for screen readers, which cannot read a canvas. */
		description: string;
	}

	let { option, height = '220px', description }: Props = $props();

	let container: HTMLDivElement;
	let chart: echarts.ECharts | null = null;

	onMount(() => {
		chart = echarts.init(container, undefined, { renderer: 'canvas' });
		chart.setOption(option, true);

		// Charts sit in a fluid layout, so track the element rather than the
		// window — a column can change width without the window doing so.
		const observer = new ResizeObserver(() => chart?.resize());
		observer.observe(container);

		return () => {
			observer.disconnect();
			chart?.dispose();
			chart = null;
		};
	});

	$effect(() => {
		// `true` replaces the option outright rather than merging, so removed
		// series do not linger from the previous render.
		chart?.setOption(option, true);
	});
</script>

<div class="chart" style:height bind:this={container} role="img" aria-label={description}></div>

<style>
	.chart {
		width: 100%;
	}
</style>
