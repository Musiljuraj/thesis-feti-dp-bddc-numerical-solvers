function save_thesis_figure(fig, filename_base)
%SAVE_THESIS_FIGURE Save figure as PDF and PNG.

  set(fig, 'paperpositionmode', 'auto');

  print(fig, [filename_base '.pdf'], '-dpdf', '-bestfit');
  print(fig, [filename_base '.png'], '-dpng', '-r300');

  fprintf('Saved: %s.pdf and %s.png\n', filename_base, filename_base);
end