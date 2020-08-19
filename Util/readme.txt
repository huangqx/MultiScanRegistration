% Go 'c++' repository
>> mex unproject.cpp linear_algebra.cpp primitive_intersect.cpp ray_mesh_intersect.cpp tri_mesh.cpp;
>> for scan_id = 1 : length(Cameras)
>>    scans{id} = simulate_scan(Shape, Cameras{scan_id});
>> end