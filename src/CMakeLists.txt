add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PUBLIC
        TOV.cpp
        TOV.hpp
        main.cpp
        constants.hpp
)

#target_compile_options(${PROJECT_NAME} PUBLIC "-O3")

