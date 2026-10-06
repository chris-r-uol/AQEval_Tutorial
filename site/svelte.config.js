import adapter from '@sveltejs/adapter-static';
import { vitePreprocess } from '@sveltejs/vite-plugin-svelte';

/** @type {import('@sveltejs/kit').Config} */
const config = {
	preprocess: vitePreprocess(),
	kit: {
		adapter: adapter(),
		paths: {
			// GitHub Pages serves a project site from /<repository>/, so the
			// deploy workflow sets BASE_PATH. Local development runs at the root.
			base: process.env.BASE_PATH ?? ''
		}
	}
};

export default config;
