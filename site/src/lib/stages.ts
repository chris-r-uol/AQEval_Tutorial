/**
 * The steps of the analysis, as the page presents them.
 *
 * A reader starts with the measurements and adds the analysis one step at a
 * time. The figure, the findings and the code all show only the steps added
 * so far. The order here is the order a reader adds them; in the code, the
 * isolation runs before the averaging (see `steps` in `code.js`).
 */
import type { Stage } from './code';

/** What the figure shows once a stage is reached. */
export const STAGE_LABEL: Record<Stage, string> = {
	data: 'Measurements',
	average: 'Averages',
	isolate: 'Isolated signal',
	breaks: 'Break points',
	changes: 'Changes'
};

/** The step that takes a reader to a stage, as an instruction. */
export const STAGE_ACTION: Record<Stage, string> = {
	data: 'Look at the measurements',
	average: 'Average the data',
	isolate: 'Isolate the signal',
	breaks: 'Find the break points',
	changes: 'Measure the changes'
};

/** What adding the step will do, in a sentence. */
export const STAGE_TEASER: Record<Stage, string> = {
	data: 'The measurements as the monitor recorded them.',
	average: 'Group the hours into blocks, so that slow changes can show through the daily ups and downs.',
	isolate: 'Take out what the weather and the seasons did, and leave what a policy could have changed.',
	breaks: 'Ask where the level of the series shifts.',
	changes: 'Fit a trend through the break points, and measure how big and how fast each change was.'
};
