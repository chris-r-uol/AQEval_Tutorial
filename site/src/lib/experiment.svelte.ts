/**
 * The state of the experiment: where every dial is set, and how far into the
 * analysis the reader has gone.
 *
 * Kept in one place because the controls, the map, the figure, the report and
 * the code panel all read the same settings. The settings are mirrored into
 * the page address, so a link reproduces an exact set-up — useful for a
 * lecturer pointing a class at one result.
 */
import { DEFAULTS, STAGES, reached, type Language, type Settings, type Stage } from './code';

const AVERAGING_CODE: Record<string, string> = { '8 hour': '8h', day: 'day', '7 day': 'week' };

/** Where a reader lands: the measurements, with no analysis yet. */
const FIRST_STAGE: Stage = 'data';

const NO_ISOLATION = { deseason: false, deweather: false, background: null } as const;

class Experiment {
	settings = $state<Settings>({ ...DEFAULTS });
	language = $state<Language>('r');
	stage = $state<Stage>(FIRST_STAGE);

	/** Change one or more dials. */
	set(change: Partial<Settings>) {
		this.settings = { ...this.settings, ...change };
	}

	/** Go to a stage: forward to add a step, back to take steps away. */
	goTo(stage: Stage) {
		this.stage = stage;
	}

	/** The stage after this one, or null at the last. */
	get nextStage(): Stage | null {
		return STAGES[STAGES.indexOf(this.stage) + 1] ?? null;
	}

	/**
	 * The settings that actually apply at this stage. Before the isolation
	 * step is added the series is analysed as measured, whatever the isolation
	 * dials are set to, so that adding the step brings the dials into play.
	 */
	get effective(): Settings {
		return reached(this.stage, 'isolate') ? this.settings : { ...this.settings, ...NO_ISOLATION };
	}

	reset() {
		this.settings = { ...DEFAULTS };
		this.stage = FIRST_STAGE;
	}

	get isDefault() {
		return (
			this.stage === FIRST_STAGE &&
			(Object.keys(DEFAULTS) as Array<keyof Settings>).every((key) => this.settings[key] === DEFAULTS[key])
		);
	}

	/** The settings and language as an address fragment. */
	toHash(): string {
		const s = this.settings;
		const params = new URLSearchParams({
			site: s.site,
			pollutant: s.pollutant,
			deseason: s.deseason ? '1' : '0',
			deweather: s.deweather ? '1' : '0',
			background: s.background ?? 'none',
			average: AVERAGING_CODE[s.averaging] ?? s.averaging,
			h: String(s.h),
			stage: this.stage,
			language: this.language
		});
		return `#${params}`;
	}

	/**
	 * Take settings from an address fragment, ignoring anything that is not a
	 * value the site actually offers: a mistyped link should fall back to the
	 * defaults rather than ask for a result that was never computed.
	 */
	fromHash(
		hash: string,
		allowed: { sites: string[]; pollutants: string[]; averaging: string[]; h: number[]; controls: string[] }
	) {
		const params = new URLSearchParams(hash.replace(/^#/, ''));
		const next: Settings = { ...DEFAULTS };

		const site = params.get('site')?.toUpperCase();
		if (site && allowed.sites.includes(site)) next.site = site;

		const pollutant = params.get('pollutant')?.toLowerCase();
		if (pollutant && allowed.pollutants.includes(pollutant)) next.pollutant = pollutant;

		if (params.has('deseason')) next.deseason = params.get('deseason') === '1';
		if (params.has('deweather')) next.deweather = params.get('deweather') === '1';

		const background = params.get('background');
		if (background === 'none') next.background = null;
		else if (background === 'air_temp') next.background = 'air_temp';
		else if (background && allowed.controls.includes(background.toUpperCase()))
			next.background = background.toUpperCase();
		// A site cannot be its own background.
		if (next.background === next.site) next.background = null;

		const average = params.get('average');
		const averaging = allowed.averaging.find((a) => (AVERAGING_CODE[a] ?? a) === average);
		if (averaging) next.averaging = averaging;

		const h = Number(params.get('h'));
		if (allowed.h.includes(h)) next.h = h;

		// If the results were rebuilt without one of the starting values, start
		// from one that exists rather than from a result that is not there.
		if (!allowed.averaging.includes(next.averaging)) next.averaging = allowed.averaging[0];
		if (!allowed.h.includes(next.h)) next.h = allowed.h[Math.floor(allowed.h.length / 2)];

		this.settings = next;
		const language = params.get('language');
		if (language === 'r' || language === 'python') this.language = language;

		// Links made before the site had stages show the whole analysis, which
		// is what they showed then. A link with no settings at all starts at
		// the beginning.
		const stage = params.get('stage') as Stage | null;
		if (stage && STAGES.includes(stage)) this.stage = stage;
		else this.stage = params.has('site') ? 'changes' : FIRST_STAGE;
	}
}

export const experiment = new Experiment();
