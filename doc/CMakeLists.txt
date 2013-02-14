#project(DECDoc)

#set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)


project(DECDoc)

#add_custom_target(doc-prerequisites
#    ALL
#    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#    )

SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )
add_custom_target(doc ALL
    COMMAND doxygen
    COMMAND ${CMAKE_COMMAND} -E rename html dec-doc
    COMMAND ${CMAKE_COMMAND} -E remove dec-doc/dec-doc.tgz
    COMMAND ${CMAKE_COMMAND} -E tar cvfz dec-doc/dec-doc.tgz dec-doc
    COMMAND ${CMAKE_COMMAND} -E rename dec-doc html
    )
