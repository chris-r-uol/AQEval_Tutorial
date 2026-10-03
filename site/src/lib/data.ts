/**
 * Loading the precomputed results.
 *
 * Nothing is calculated in the browser. `pipeline/build_data.py` ran the real
 * analysis for every combination of dials and wrote the answers to
 * `static/data`; moving a dial fetches the file that holds the answer. Files
 * are small and are kept once fetched, so going back to a setting is instant.
 */
import { base } from '$app/paths';
import type { Language, Settings } from './code';
import type { Manifest, RunFile } from './types';
import type { Feature, Polygon } from 'geojson';

const cache = new Map<string, Promise<unknown>>();

function load<T>(path: string): Promise<T> {
	let pending = cache.get(path);
	if (!pending) {
		pending = fetch(`${base}/data/${path}`).then((response) => {
			if (!response.ok) throw new Error(`Could not load ${path} (${response.status})`);
			return response.json();
		});
		// A failed request should be retried next time, not remembered.
		pending.catch(() => cache.delete(path));
		cache.set(path, pending);
	}
	return pending as Promise<T>;
}

export const loadManifest = () => load<Manifest>('manifest.json');
export const loadBoundary = () => load<Feature<Polygon>>('caz.geojson');

/** File-name key for one isolation setting. Must match `iso_key` in the pipeline. */
export function isolationKey(settings: Pick<Settings, 'deseason' | 'deweather' | 'background'>) {
	const background = (settings.background ?? 'none').toLowerCase();
	return `ds${Number(settings.deseason)}-dw${Number(settings.deweather)}-bg${background}`;
}

/**
 * Every result for one site, pollutant and isolation setting.
 *
 * R and Python do not give identical answers, so each language has its own
 * complete set of results, computed by that language.
 */
export const loadRuns = (settings: Settings, language: Language) =>
	load<RunFile>(`runs/${language}/${settings.site}_${settings.pollutant}_${isolationKey(settings)}.json`);

/**
 * The results with nothing isolated. Its series is the measurements
 * themselves, which the figure draws underneath an isolated signal.
 */
export const loadMeasured = (settings: Settings, language: Language) =>
	loadRuns({ ...settings, deseason: false, deweather: false, background: null }, language);
