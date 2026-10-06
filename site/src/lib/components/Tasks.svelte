<script lang="ts">
	/**
	 * The practical's tasks.
	 *
	 * Each can be answered on this page, by moving the dials. Each also names
	 * what to change in the code, for doing it for real in the workspace.
	 */
	import { CODESPACE_URL } from '$lib/config';
</script>

<ol class="tasks">
	<li>
		<h3>Does NOx tell the same story as NO₂?</h3>
		<p>
			Measure the changes in NOx at Bradford Mayo Avenue. Pin the NO₂ result first, then switch
			pollutant. Do the changes fall on the same dates? Are they the same size?
		</p>
		<p class="muted">In the code: change <code>"no2"</code> to <code>"nox"</code>.</p>
	</li>
	<li>
		<h3>Did it happen where there was no Clean Air Zone?</h3>
		<p>
			Repeat the analysis for NO₂ and NOx at Leeds Headingley Kerbside, a roadside site with no charging
			zone. If a change shows up at both sites on the same date, what does that say about its cause?
		</p>
		<p class="muted">In the code: change the site from <code>"bdma"</code> to <code>"led6"</code>.</p>
	</li>
	<li>
		<h3>Which site makes the better control?</h3>
		<p>
			Leeds Centre and Dewsbury Ashworth Grange are both background sites. Use each in turn under
			<em>Also account for</em>. Which changes survive, and which disappear? Look at the site comparison
			chart too: how complete is each record, and how far away is each site?
		</p>
		<p class="muted">
			In the code: a second download, joined to the first by date, then passed as <code>background</code>.
		</p>
	</li>
	<li>
		<h3>How much does the sensitivity matter?</h3>
		<p>
			For NO₂ at Bradford Mayo Avenue, make the analysis more sensitive: move <code>h</code> from 0.3 down
			to 0.12, pinning as you go. Which changes appear at every setting, and which come and go? Which
			result would you report, and how would you justify the choice?
		</p>
		<p class="muted">In the code: change <code>h</code>. Smaller values take longer to run.</p>
	</li>
	<li>
		<h3>Take it further</h3>
		<p>
			The page only knows the set-ups it was given. In the <a href={CODESPACE_URL} target="_blank" rel="noopener"
				>workspace</a
			> you can run anything: another site, other years, another pollutant, or your own choice of background.
			Pick a UK city with a Clean Air Zone and find out whether its air changed.
		</p>
	</li>
</ol>

<style>
	.tasks {
		margin: 0;
		padding: 0;
		list-style: none;
		counter-reset: task;
		display: grid;
		grid-template-columns: repeat(auto-fit, minmax(min(100%, 360px), 1fr));
		gap: var(--space-3);
	}

	.tasks li {
		counter-increment: task;
		padding: var(--space-4);
		background: var(--surface);
		border-radius: var(--radius-lg);
		box-shadow: var(--shadow-raised);
		border-top: 4px solid var(--colour-brand-dark);
	}

	.tasks h3 {
		font-size: 1rem;
		margin-bottom: var(--space-2);
	}

	.tasks h3::before {
		content: 'Task ' counter(task) '. ';
		color: var(--text-muted);
	}

	.tasks p {
		font-size: 0.9375rem;
	}

	.tasks p:last-child {
		margin-bottom: 0;
	}
</style>
