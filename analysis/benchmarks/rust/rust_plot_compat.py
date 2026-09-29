"""Original upstream plotting function with one Python-3 iterator compatibility fix.
No RUST numerical calculation is replaced or changed.
"""
import inspect


def compatible_metafootprint(module):
    source = inspect.getsource(module.RUST_metagene_plot)
    old = 'coverage = map(float, linesplit[2:])'
    new = 'coverage = list(map(float, linesplit[2:]))'
    if old not in source and new in source:
        return module.RUST_metagene_plot
    if source.count(old) != 1:
        raise RuntimeError('Upstream plot function changed; inspect compatibility adapter')
    namespace = dict(module.__dict__)
    exec(compile(source.replace(old, new), '<RUST plot iterator compatibility>', 'exec'), namespace)
    return namespace['RUST_metagene_plot']


def render_original_plot(module, profile, output):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.69, 6.0))
    with profile.open() as handle:
        compatible_metafootprint(module)(handle, ax)
    fig.savefig(output, dpi=180, bbox_inches='tight')
    plt.close(fig)
