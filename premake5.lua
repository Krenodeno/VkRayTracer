-- premake5.lua
workspace "compute"

	configurations { "Debug", "DebugSAN", "Release" }
	architecture "x64"

	location "build/"
	debugdir "./"
	targetdir "bin/"
	objdir "build/obj/"

	local gKitDir = "/home/thomas/Libraries/gkit2light/"
	local gKitFiles = { gKitDir .. "src/gKit/*.cpp",  gKitDir .. "src/gKit/*.h" }

	includedirs (gKitDir .. "src/gKit")

	links { "SDL2", "SDL2_image", "GLEW", "GL" }

	language "C++"
	cppdialect "C++11"

	includedirs { "extern/include/" }

	filter "configurations:Debug*"
		symbols "On"
		optimize "Debug"

	filter "configurations:Release"
		symbols "Off"
		optimize "On"
		buildoptions { "-fopenmp" }
		linkoptions { "-fopenmp" }

	filter "configurations:*SAN"
		buildoptions { "-fno-omit-frame-pointer" }
		buildoptions { "-fsanitize=undefined,address" }
		linkoptions { "-fsanitize=undefined,address" }

	filter {}

project "gKit"
	kind "StaticLib"
	files (gKitFiles)

project "GPURayTracing"
	kind "ConsoleApp"
	files { "src/ComputeApp.hpp", "src/ComputeApp.cpp", "src/computeRayTracing.cpp" }

	links { "gKit", "vulkan" }

project "ray_tuto"
	kind "ConsoleApp"
	files { "src/ray_tuto.cpp" }

	links "gKit"

project "ray_tuto_mc"
	kind "ConsoleApp"
	files { "src/ray_tuto_mc.cpp" }

	links "gKit"

project "tuto_bvh"
	kind "ConsoleApp"
	files { "src/tuto_bvh.cpp" }

	links "gKit"

project "tuto_bvh_simple"
	kind "ConsoleApp"
	files { "src/tuto_bvh_simple.cpp" }

	links "gKit"
