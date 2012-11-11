#============================================================================
#
# Program: Bender
#
# Copyright (c) Kitware Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#============================================================================

cmake_minimum_required(VERSION 2.8.2)

set(kit Bender)
project(${kit})

#-----------------------------------------------------------------------------
# Library mode: SHARED (default) or STATIC
option(BUILD_SHARED_LIBS "Build Bender with shared libraries." ON)
mark_as_advanced(BUILD_SHARED_LIBS)
set(${kit}_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0016)
  cmake_policy(SET CMP0016 NEW)
endif()
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#
option(${kit}_SUPERBUILD "Build ${kit} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${kit}_SUPERBUILD)

#-----------------------------------------------------------------------------
# Project version number.
#
set(${kit}_VERSION_MAJOR 0)
set(${kit}_VERSION_MINOR 1)
set(${kit}_VERSION_PATCH 0)
set(${kit}_VERSION
    "${${kit}_VERSION_MAJOR}.${${kit}_VERSION_MINOR}.${${kit}_VERSION_PATCH}")

#-----------------------------------------------------------------------------
# Update CMake module path
# Note: FindXXX.cmake script specific to utility should be copied into Utilities/CMake
#
set(CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake"
  ${CMAKE_MODULE_PATH})

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
find_package(Git)
if(NOT GIT_FOUND)
  message(FATAL_ERROR "error: Install Git and try to re-configure")
endif()

#-----------------------------------------------------------------------------
# Superbuild script
#

if(${kit}_SUPERBUILD)
  add_subdirectory("CMake/SuperBuild")
  return()
endif()

