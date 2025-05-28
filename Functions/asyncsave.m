function asyncsave(root, filename, struct)
    while exist(fullfile(root, '.lock'), 'file') > 0
        pause(0.1);
    end
    save(fullfile(root, '.lock'));
    save(fullfile(root, filename)', '-struct', 'struct', '-mat');
    delete(fullfile(root, '.lock'));
end