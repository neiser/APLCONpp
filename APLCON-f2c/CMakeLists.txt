INCLUDE_DIRECTORIES(.)
add_library(aplcon-f2c aploop.F.c)
target_link_libraries(aplcon-f2c f2c m)