"""Posterior mean/median tables with asymmetric equal-tail interval errors."""
from html import escape
import numpy as np
from sample_plot_labels import plot_label
from posterior_display import interval_text, number


def posterior_statistics_html(data, source):
    samples=np.asarray(data['samples'],dtype=float)
    labels=data['labels']
    options=data['options']
    rounded=options.get('round_statistics',True)
    show_mad=options.get('show_mad',False)
    level=float(options.get('quantiles') or 68.3)
    tail=(100.-level)/2
    rows=[f'<h1>{escape(source)} posterior statistics</h1>',
          f'<p>{len(samples)} posterior rows; {len(labels)} parameters. '
          f'Errors give the distances to the {level:g}% equal-tail credible interval '
          '(not standard errors of the mean). Periodic angles use a continuous '
          'branch around their circular mean. Finite values are used per parameter.</p>']
    populations=[('All posterior samples',np.ones(len(samples),dtype=bool))]
    if data.get('stable') is not None:
        populations.append(('Stable posterior samples',np.asarray(data['stable'],dtype=bool)))
    for title,mask in populations:
        rows.append(f'<h2>{title} ({int(mask.sum())} rows)</h2>')
        if not mask.any():
            rows.append('<p>No samples satisfy the stability threshold.</p>')
            continue
        headings=['Parameter','Finite rows','Mean − / +','Median − / +','Standard deviation']
        if show_mad:headings.append('MAD')
        rows.append('<table border="1" cellspacing="0" cellpadding="6"><tr>'+''.join(f'<th>{h}</th>' for h in headings)+'</tr>')
        for j,label in enumerate(labels):
            x=samples[mask,j];x=x[np.isfinite(x)]
            cells=[plot_label(label),str(len(x))]
            if len(x):
                low,high=np.percentile(x,[tail,100-tail]);scale=(high-low)/2
                cells += [interval_text(np.mean(x),low,high,rounded),
                          interval_text(np.median(x),low,high,rounded),
                          number(np.std(x,ddof=1),scale,rounded) if len(x)>1 else 'undefined']
                if show_mad:cells.append(number(np.median(np.abs(x-np.median(x))),scale,rounded))
            else:cells+=['undefined']*(len(headings)-2)
            rows.append('<tr>'+''.join(f'<td>{c}</td>' for c in cells)+'</tr>')
        rows.append('</table>')
    return ''.join(rows)
