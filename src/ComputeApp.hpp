#ifndef COMPUTE_APP_HPP
#define COMPUTE_APP_HPP

#include <vulkan/vulkan.hpp>

#include <cstdio>
#include <cmath>

#include "image.h"

/** Callback */

static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
	VkDebugUtilsMessageSeverityFlagBitsEXT,
	VkDebugUtilsMessageTypeFlagsEXT,
	const VkDebugUtilsMessengerCallbackDataEXT*,
	void*);

static std::vector<uint32_t> readShaderFile(uint32_t& length, std::string filename) {
	FILE* fp = fopen(filename.c_str(), "rb");
	if (fp == NULL) {
		throw std::runtime_error("Could not find or open file\n");
	}

	// get file size.
	fseek(fp, 0, SEEK_END);
	long filesize = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	long filesizepadded = long(ceil(filesize / 4.0)) * 4;

	// read file contents.
	std::vector<uint32_t> str(filesizepadded);
	fread(str.data(), filesize, sizeof(char), fp);
	fclose(fp);

	// data padding.
	for (int i = filesize; i < filesizepadded; i++) {
		str[i] = 0;
	}

	length = filesizepadded;
	return str;
}

const int WORKGROUP_SIZE = 32; // Workgroup size in compute shader.

class ComputeApp {
public:
	ComputeApp();
	ComputeApp(int im_width, int im_height);

	~ComputeApp();

	void run();

	void addBuffer(uint64_t bufferSize);
	void fillBuffer(uint32_t index, const void* data, size_t dataSize);

	template<typename T>
	std::vector<T> getDataFromBuffer(uint32_t index, size_t dataSize) {
		// Get data back from Device Memory
		auto data = static_cast<T*>(device.mapMemory(buffersMemory[index], /*offset*/ 0, dataSize));
		std::vector<T> res(data, data + dataSize);
		device.unmapMemory(buffersMemory[index]);
		return res;
	}

private:

	int width, height;

	struct QueueFamilyIndices;

	vk::Instance instance;
	vk::DebugUtilsMessengerEXT callback;

	vk::PhysicalDevice physicalDevice;

	vk::Device device;
	vk::Queue computeQueue;

	vk::Pipeline pipeline;
	vk::PipelineLayout pipelineLayout;
	vk::ShaderModule computeShaderModule;

	vk::CommandPool commandPool;
	vk::CommandBuffer commandBuffer;

	vk::DescriptorPool descriptorPool;
	std::vector<vk::DescriptorSet> descriptorSets;
	vk::DescriptorSetLayout descriptorSetLayout;

	std::vector<vk::Buffer> buffers;
	std::vector<vk::DeviceMemory> buffersMemory;
	std::vector<uint64_t> bufferSizes;

	void init();

	void draw();

	void cleanup();

	void createInstance();
	void setupDebugCallback();
	void pickPhysicalDevice();
	void createDevice();
	void createBuffer(vk::DeviceSize, vk::BufferUsageFlags, vk::MemoryPropertyFlags, vk::Buffer&, vk::DeviceMemory&);
	void createBuffers();
	void createDescriptorSetLayout();
	void createDescriptorSet();
	void createComputePipeline();
	void createCommandeBuffer();

	bool isDeviceSuitable(vk::PhysicalDevice);

	QueueFamilyIndices findQueueFamilies(vk::PhysicalDevice);

	uint32_t findMemoryType(uint32_t, vk::MemoryPropertyFlags);

	struct QueueFamilyIndices {
		int computeFamily = -1;

		bool isComplete() {
			return computeFamily >= 0;
		}
	};

	struct UniformBufferObject {
		unsigned int width;
		unsigned int height;
	};
};

#endif
