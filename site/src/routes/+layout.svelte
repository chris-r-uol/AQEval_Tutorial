<script lang="ts">
	import '$lib/theme.css';
	import { base } from '$app/paths';
	import { CODESPACE_URL, REPOSITORY_URL } from '$lib/config';

	let { children } = $props();
</script>

<svelte:head>
	<link rel="icon" href="{base}/favicon.svg" />
	<title>AQEval tutorial: finding changes in Bradford's air</title>
	<meta
		name="description"
		content="An interactive tutorial on detecting changes in air quality with AQEval, using monitoring data from Bradford."
	/>
</svelte:head>

<a class="skip" href="#main">Skip to the experiment</a>

<header class="top">
	<a class="top__brand" href="{base}/">AQEval tutorial</a>
	<nav class="top__nav" aria-label="Sections">
		<a href="#experiment">Experiment</a>
		<a href="#tasks">Tasks</a>
		<a href="#about">About</a>
	</nav>
	<div class="top__links">
		<a href={REPOSITORY_URL} target="_blank" rel="noopener">GitHub</a>
		<a class="top__run" href={CODESPACE_URL} target="_blank" rel="noopener">Run the code</a>
	</div>
</header>

<main id="main">
	{@render children()}
</main>

<style>
	.skip {
		position: absolute;
		left: -9999px;
		z-index: 2000;
		padding: var(--space-2) var(--space-4);
		background: var(--focus);
		color: var(--text);
	}

	.skip:focus {
		left: 0;
		top: 0;
	}

	.top {
		position: sticky;
		top: 0;
		z-index: 1000;
		display: flex;
		align-items: center;
		gap: var(--space-5);
		height: var(--header-height);
		padding: 0 var(--space-4);
		background: var(--colour-brand-dark);
		color: var(--text-inverse);
	}

	.top a {
		color: var(--text-inverse);
		text-decoration: none;
	}

	.top a:hover {
		text-decoration: underline;
	}

	.top__brand {
		font-weight: 700;
		font-size: 1.125rem;
		white-space: nowrap;
	}

	.top__nav {
		display: flex;
		gap: var(--space-4);
		flex: 1;
		font-size: 0.9375rem;
	}

	.top__links {
		display: flex;
		align-items: center;
		gap: var(--space-4);
		font-size: 0.9375rem;
	}

	.top__run {
		padding: var(--space-1) var(--space-3);
		background: var(--surface);
		color: var(--colour-brand-dark) !important;
		border-radius: var(--radius);
		font-weight: 700;
		white-space: nowrap;
	}

	@media (max-width: 640px) {
		.top__nav {
			display: none;
		}

		.top {
			justify-content: space-between;
		}
	}
</style>
