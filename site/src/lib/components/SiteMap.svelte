<script lang="ts">
	/**
	 * The map: Bradford, its Clean Air Zone, and the monitoring sites.
	 *
	 * Leaflet with raster OpenStreetMap tiles — no API key and no tile account.
	 * The sites are labelled with a figure for the chosen pollutant, so the map
	 * is a comparison in its own right: which sites are high, and which fell.
	 */
	import { onMount } from 'svelte';
	// Type-only, so it is erased at build time. The runtime import happens in
	// onMount: Leaflet touches `window` at module scope.
	import type * as Leaflet from 'leaflet';
	import 'leaflet/dist/leaflet.css';
	import type { Feature, Polygon } from 'geojson';
	import { SITE_COLOURS } from '$lib/charts';
	import { UNIT, formatNumber, formatSigned } from '$lib/format';
	import type { MapMetric, Site } from '$lib/types';

	interface Props {
		sites: Site[];
		boundary: Feature<Polygon> | null;
		pollutant: string;
		selected: string;
		/** The site being used as the background term, if one is. */
		background: string | null;
		metric: MapMetric;
		onSelect: (code: string) => void;
	}

	let { sites, boundary, pollutant, selected, background, metric, onSelect }: Props = $props();

	let container: HTMLDivElement;
	// Raw: Svelte must not wrap the Leaflet library in a reactive proxy.
	let L = $state.raw<typeof Leaflet | null>(null);
	let map: Leaflet.Map | null = null;
	let siteLayer: Leaflet.LayerGroup | null = null;
	let boundaryLayer: Leaflet.LayerGroup | null = null;

	onMount(() => {
		let disposed = false;

		// Async work inside a synchronous onMount, so the teardown is registered
		// with Svelte immediately rather than after the import resolves.
		(async () => {
			const leaflet = await import('leaflet');
			if (disposed) return;
			const library = (leaflet.default ?? leaflet) as typeof Leaflet;

			map = library.map(container, {
				// A map that zooms as the page scrolls past it traps the reader.
				scrollWheelZoom: false,
				zoomControl: true,
				attributionControl: true
			});
			map.attributionControl.setPrefix('');
			library
				.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
					maxZoom: 19,
					attribution: '© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
				})
				.addTo(map);
			library.control.scale({ metric: true, imperial: false }).addTo(map);

			boundaryLayer = library.layerGroup().addTo(map);
			siteLayer = library.layerGroup().addTo(map);

			narrow = container.clientWidth < 520;
			// Room at the sides for the labels, which sit beside the markers.
			map.fitBounds(library.latLngBounds(sites.map((site) => [site.lat, site.lon])), {
				padding: [narrow ? 64 : 110, 44]
			});
			// Setting L last triggers the effects below, now the map exists.
			L = library;
		})();

		// The map sits in a fluid column; Leaflet only learns of a resize if told.
		const observer = new ResizeObserver(() => {
			map?.invalidateSize();
			narrow = container.clientWidth < 520;
		});
		observer.observe(container);

		return () => {
			disposed = true;
			observer.disconnect();
			map?.remove();
			map = null;
		};
	});

	function value(site: Site): number | null {
		const stats = site.stats[pollutant];
		if (metric === 'before') return stats.before;
		if (metric === 'after') return stats.after;
		if (stats.before === null || stats.after === null || stats.before === 0) return null;
		return ((stats.after - stats.before) / stats.before) * 100;
	}

	// On a phone the full names run off the edge of the map, so the labels fall
	// back to the site codes, which the list of sites also shows.
	let narrow = $state(false);

	function label(site: Site): string {
		const figure = value(site);
		const text = metric === 'change' ? formatSigned(figure, 0, '%') : `${formatNumber(figure)} ${UNIT}`;
		return `<strong>${narrow ? site.code : site.name}</strong><br>${text}`;
	}

	$effect(() => {
		if (!L || !boundaryLayer) return;
		boundaryLayer.clearLayers();
		if (boundary) {
			L.geoJSON(boundary, {
				interactive: false,
				style: { color: '#1d70b8', weight: 2, fillColor: '#1d70b8', fillOpacity: 0.1 }
			}).addTo(boundaryLayer);
		}
	});

	$effect(() => {
		if (!L || !siteLayer) return;
		const library = L;
		siteLayer.clearLayers();

		// A dashed line ties the site to the one it is being compared against.
		const from = sites.find((site) => site.code === selected);
		const to = sites.find((site) => site.code === background);
		if (from && to) {
			library
				.polyline(
					[
						[from.lat, from.lon],
						[to.lat, to.lon]
					],
					{ color: '#0b0c0c', weight: 1.5, dashArray: '5 5', interactive: false }
				)
				.addTo(siteLayer);
		}

		// Labels alternate sides so neighbouring sites do not print over each other.
		const west = [...sites].sort((a, b) => a.lon - b.lon);
		for (const site of sites) {
			const isSelected = site.code === selected;
			const marker = library
				.circleMarker([site.lat, site.lon], {
					radius: isSelected ? 11 : 8,
					color: isSelected ? '#0b0c0c' : '#ffffff',
					weight: isSelected ? 3 : 2,
					fillColor: SITE_COLOURS[site.type] ?? '#505a5f',
					fillOpacity: 1
				})
				.bindTooltip(label(site), {
					permanent: true,
					direction: west.indexOf(site) % 2 === 0 ? 'left' : 'right',
					offset: [west.indexOf(site) % 2 === 0 ? -10 : 10, 0],
					className: isSelected ? 'site-label site-label--on' : 'site-label'
				})
				.on('click', () => onSelect(site.code))
				.addTo(siteLayer!);
			// Markers are focusable SVG paths; say what they are and let the
			// keyboard choose them too.
			const element = marker.getElement();
			if (element) {
				element.setAttribute('tabindex', '0');
				element.setAttribute('role', 'button');
				element.setAttribute('aria-label', `Analyse ${site.name}`);
				element.addEventListener('keydown', (event) => {
					const key = (event as KeyboardEvent).key;
					if (key === 'Enter' || key === ' ') {
						event.preventDefault();
						onSelect(site.code);
					}
				});
			}
		}
	});
</script>

<div class="map">
	<div class="map__canvas" bind:this={container}></div>
	<ul class="map__legend" aria-label="Map key">
		<li><span class="dot" style:background={SITE_COLOURS['Urban Traffic']}></span> Roadside site</li>
		<li><span class="dot" style:background={SITE_COLOURS['Urban Background']}></span> Background site</li>
		<li><span class="zone"></span> Clean Air Zone</li>
	</ul>
</div>

<style>
	.map {
		position: relative;
		border-radius: var(--radius);
		overflow: hidden;
		/* Keeps Leaflet's stacking inside the card, below the sticky header. */
		isolation: isolate;
	}

	.map__canvas {
		height: 340px;
		background: var(--surface-sunken);
	}

	.map__legend {
		position: absolute;
		left: var(--space-2);
		bottom: var(--space-5);
		z-index: 500;
		margin: 0;
		padding: var(--space-2) var(--space-3);
		list-style: none;
		background: rgba(255, 255, 255, 0.94);
		border-radius: var(--radius);
		box-shadow: var(--shadow-raised);
		font-size: 0.75rem;
		line-height: 1.6;
	}

	.dot {
		display: inline-block;
		width: 0.625rem;
		height: 0.625rem;
		border-radius: 50%;
		margin-right: 2px;
	}

	.zone {
		display: inline-block;
		width: 0.875rem;
		height: 0.625rem;
		border: 2px solid var(--colour-brand);
		background: rgba(29, 112, 184, 0.15);
		margin-right: 2px;
		vertical-align: -1px;
	}

	:global(.leaflet-container) {
		font: inherit;
	}

	:global(.leaflet-control-attribution) {
		font-size: 0.6875rem;
	}

	:global(.leaflet-interactive:focus-visible) {
		outline: 3px solid var(--focus);
	}

	/* Created by Leaflet, outside this component's scoped styles. */
	:global(.site-label) {
		font-size: 0.75rem;
		line-height: 1.3;
		padding: 2px var(--space-2);
		border: 1px solid var(--border);
		border-radius: var(--radius);
		box-shadow: none;
		color: var(--text);
	}

	:global(.site-label--on) {
		border: 2px solid var(--text);
	}

	:global(.site-label::before) {
		display: none;
	}
</style>
