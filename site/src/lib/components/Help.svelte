<script lang="ts">
	/** A "what does this do?" disclosure under a dial. */
	import type { Snippet } from 'svelte';

	interface Props {
		summary?: string;
		children: Snippet;
	}

	let { summary = 'What does this do?', children }: Props = $props();
</script>

<details class="help">
	<summary>{summary}</summary>
	<div class="help__body">{@render children()}</div>
</details>

<style>
	.help {
		margin-top: var(--space-2);
		font-size: 0.875rem;
	}

	summary {
		color: var(--colour-brand);
		cursor: pointer;
		text-decoration: underline;
		width: fit-content;
	}

	summary:hover {
		color: var(--colour-brand-dark);
	}

	.help__body {
		margin-top: var(--space-2);
		padding-left: var(--space-3);
		border-left: 4px solid var(--border-subtle);
		color: var(--text);
	}

	.help__body :global(p:last-child) {
		margin-bottom: 0;
	}
</style>
