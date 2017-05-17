function filename_ell = empty2nan(filename_ell)      %named function
  filename_ell(cellfun(@isempty, filename_ell)) = ['nan']
end