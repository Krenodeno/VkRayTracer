-- premake5.lua
workspace "compute"

	configurations { "Debug", "Release" }
	platforms "x64"
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

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"
		buildoptions { "-fopenmp" }
		linkoptions { "-fopenmp" }

	filter {}

project "gKit"
	kind "StaticLib"
	files (gKitFiles)

project "ray_tuto"
	kind "ConsoleApp"
	files { "src/ray_tuto.cpp" }

	links "gKit"

project "ray_tuto_mc"
	kind "ConsoleApp"
	files { "src/ray_tuto_mc.cpp" }

	links "gKit"
