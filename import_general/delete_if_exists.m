function delete_if_exists(file)
    if exist(file, 'file') == 2
        delete(file)
    end
end