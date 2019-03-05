#include "ComputeApp.hpp"

#include <cassert>
#include <fstream>
#include <iostream>

static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
	VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
	VkDebugUtilsMessageTypeFlagsEXT messageType,
	const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData,
	void* pUserData) {

	std::ostream& out = [&]() -> std::ostream& {
		if (messageSeverity == VkDebugUtilsMessageSeverityFlagBitsEXT::VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT)
			return (std::cout);
		if (messageSeverity == VkDebugUtilsMessageSeverityFlagBitsEXT::VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT)
			return (std::cerr);
		return (std::clog);
	}();

	if (messageType & VkDebugUtilsMessageTypeFlagBitsEXT::VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT)
		out << "[GENERAL]     ";
	if (messageType & VkDebugUtilsMessageTypeFlagBitsEXT::VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT)
		out << "[PERFORMANCE] ";
	if (messageType & VkDebugUtilsMessageTypeFlagBitsEXT::VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT)
		out << "[VALIDATION]  ";

	out << pCallbackData->pMessage << "\n";

	return VK_FALSE;
}

ComputeApp::ComputeApp() {

}

ComputeApp::~ComputeApp() {
	cleanup();
}

void ComputeApp::run() {
	init();

	draw();
}

void ComputeApp::draw() {
	vk::SubmitInfo submitInfo;
	submitInfo.commandBufferCount = 1;
	submitInfo.pCommandBuffers = &commandBuffer;

	vk::Fence fence;
	vk::FenceCreateInfo fenceCreateInfo;
	fence = device.createFence(fenceCreateInfo);

	computeQueue.submit(submitInfo, fence);

	device.waitForFences(fence, VK_TRUE, std::numeric_limits<uint64_t>::max());

	device.destroyFence(fence);
}

int ComputeApp::addBuffer(uint64_t bufferSize) {
	bufferSizes.push_back(bufferSize);
	return (int)bufferSizes.size() - 1;
}

void ComputeApp::fillBuffer(uint32_t index, const void* dataToCopy, uint64_t dataSize) {
	auto data = device.mapMemory(buffersMemory[index], /*offset*/ 0, dataSize, vk::MemoryMapFlags());
		std::memcpy(data, dataToCopy, static_cast<size_t>(dataSize));
	device.unmapMemory(buffersMemory[index]);
}


void ComputeApp::init() {
	createInstance();
	setupDebugCallback();
	pickPhysicalDevice();
	createDevice();
	createBuffers();
	createDescriptorSetLayout();
	createDescriptorSet();
	createComputePipeline();
	createCommandeBuffer();
}

void ComputeApp::cleanup() {
	// Destroy vulkan objects
	device.destroyShaderModule(computeShaderModule);
	device.destroyCommandPool(commandPool);
	device.destroyPipeline(pipeline);
	device.destroyPipelineLayout(pipelineLayout);
	device.destroyDescriptorPool(descriptorPool);
	device.destroyDescriptorSetLayout(descriptorSetLayout);

	for (int i = 0; i < buffers.size(); i ++) {
		device.freeMemory(buffersMemory[i]);
		device.destroyBuffer(buffers[i]);
	}
	// Then destroy Device
	device.destroy();
	// destroy debug utils
	auto func = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>(instance.getProcAddr("vkDestroyDebugUtilsMessengerEXT"));
	func(instance, static_cast<VkDebugUtilsMessengerEXT>(callback), nullptr);
	// Finally destroy the instance
	instance.destroy();
}

void ComputeApp::createInstance() {
	vk::ApplicationInfo appInfo;
	appInfo.pApplicationName = "Compute Application";
	appInfo.applicationVersion = VK_MAKE_VERSION(0, 0, 0);
	appInfo.pEngineName = "NOENGINE";
	appInfo.engineVersion = VK_MAKE_VERSION(0, 0, 0);
	appInfo.apiVersion = VK_API_VERSION_1_0;

	std::vector<const char*> extensions = { VK_EXT_DEBUG_UTILS_EXTENSION_NAME };
	std::vector<const char*> layers = { "VK_LAYER_LUNARG_standard_validation" };

	vk::InstanceCreateInfo createInfo;
	createInfo.pApplicationInfo = &appInfo;
	createInfo.enabledExtensionCount = extensions.size();
	createInfo.ppEnabledExtensionNames = extensions.data();
	createInfo.enabledLayerCount = layers.size();
	createInfo.ppEnabledLayerNames = layers.data();

	instance = vk::createInstance(createInfo);
}

void ComputeApp::setupDebugCallback() {
	vk::DebugUtilsMessengerCreateInfoEXT createInfo;
	createInfo.messageSeverity =
		//vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
		//vk::DebugUtilsMessageSeverityFlagBitsEXT::eInfo |
		vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
		vk::DebugUtilsMessageSeverityFlagBitsEXT::eError;
	createInfo.messageType =
		//vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
		vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
		vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance;
	createInfo.pfnUserCallback = debugCallback;
	createInfo.pUserData = nullptr;

	auto func = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(instance.getProcAddr("vkCreateDebugUtilsMessengerEXT"));
	VkDebugUtilsMessengerEXT tmp;
	auto info = (VkDebugUtilsMessengerCreateInfoEXT)createInfo;
	func(instance, &info, nullptr, &tmp);
	callback = tmp;
}

void ComputeApp::pickPhysicalDevice() {

	auto devices = instance.enumeratePhysicalDevices();

	if (devices.empty())
		throw std::runtime_error("Failed to find GPUs with Vulkan suport!");

	for (const auto& device: devices) {
		if (isDeviceSuitable(device)) {
			physicalDevice = device;
			break;
		}
	}

	auto properties = physicalDevice.getProperties();
	auto deviceLimits = properties.limits;

	std::cout << "Limits : Max Compute Workgroup Count : x=" << deviceLimits.maxComputeWorkGroupCount[0];
	std::cout << ", y=" << deviceLimits.maxComputeWorkGroupCount[1];
	std::cout << ", z=" << deviceLimits.maxComputeWorkGroupCount[2];
	std::cout << "\nMax Compute Workgroup Invocations : " << deviceLimits.maxComputeWorkGroupInvocations;
	std::cout << "\nMax Compute Workgroup Size : x=" << deviceLimits.maxComputeWorkGroupSize[0];
	std::cout << ", y=" << deviceLimits.maxComputeWorkGroupSize[1];
	std::cout << ", z=" << deviceLimits.maxComputeWorkGroupSize[2];
	std::cout << std::endl;

	if (!physicalDevice)
		throw std::runtime_error("Failed to find a suitable GPU!");
}

void ComputeApp::createDevice() {
	QueueFamilyIndices indices = findQueueFamilies(physicalDevice);

	// Queues
	float priority = 1.f;
	vk::DeviceQueueCreateInfo queueCreateInfo;
	queueCreateInfo.queueFamilyIndex = indices.computeFamily;
	queueCreateInfo.queueCount = 1;
	queueCreateInfo.pQueuePriorities = &priority;


	vk::PhysicalDeviceFeatures physicalDeviceFeatures;

	std::vector<const char*> extensions = {};

	vk::DeviceCreateInfo createInfo;
	createInfo.queueCreateInfoCount = 1;
	createInfo.pQueueCreateInfos = &queueCreateInfo;
	createInfo.pEnabledFeatures = &physicalDeviceFeatures;
	createInfo.enabledExtensionCount = extensions.size();
	createInfo.ppEnabledExtensionNames = extensions.data();
	createInfo.enabledLayerCount = 0;

	device = physicalDevice.createDevice(createInfo);

	// get Queue
	computeQueue = device.getQueue(indices.computeFamily, 0);

}

void ComputeApp::createBuffer(vk::DeviceSize size, vk::BufferUsageFlags usage, vk::MemoryPropertyFlags properties, vk::Buffer& buffer, vk::DeviceMemory& memory) {
	vk::BufferCreateInfo bufferInfo;
	bufferInfo.size = size;
	bufferInfo.usage = usage;
	bufferInfo.sharingMode = vk::SharingMode::eExclusive;

	buffer = device.createBuffer(bufferInfo);

	vk::MemoryRequirements memRequirements = device.getBufferMemoryRequirements(buffer);

	vk::MemoryAllocateInfo allocInfo;
	allocInfo.allocationSize = memRequirements.size;
	allocInfo.memoryTypeIndex = findMemoryType(memRequirements.memoryTypeBits, properties);

	memory = device.allocateMemory(allocInfo);

	device.bindBufferMemory(buffer, memory, /*MemoryOffset*/ 0);

}

void ComputeApp::createBuffers() {
	for (int i = 0; i < bufferSizes.size(); i++) {
		buffers.push_back(vk::Buffer());
		buffersMemory.push_back(vk::DeviceMemory());
		createBuffer(bufferSizes[i], vk::BufferUsageFlagBits::eStorageBuffer, vk::MemoryPropertyFlagBits::eHostCoherent | vk::MemoryPropertyFlagBits::eHostVisible, buffers[i], buffersMemory[i]);
	}
}

void ComputeApp::createDescriptorSetLayout() {
	std::vector<vk::DescriptorSetLayoutBinding> computeBindings(buffers.size());

	for (int i = 0; i < computeBindings.size(); i++) {
		computeBindings[i].binding = i;
		computeBindings[i].descriptorType = vk::DescriptorType::eStorageBuffer;
		computeBindings[i].descriptorCount = 1;
		computeBindings[i].stageFlags = vk::ShaderStageFlagBits::eCompute;
	}

	vk::DescriptorSetLayoutCreateInfo createInfo;
	createInfo.bindingCount = computeBindings.size();
	createInfo.pBindings = computeBindings.data();

	descriptorSetLayout = device.createDescriptorSetLayout(createInfo);
}

void ComputeApp::createDescriptorSet() {
	vk::DescriptorPoolSize poolSize;
	poolSize.type = vk::DescriptorType::eStorageBuffer;
	poolSize.descriptorCount = buffers.size();

	vk::DescriptorPoolCreateInfo createInfo;
	createInfo.maxSets = 2;
	createInfo.poolSizeCount = 1;
	createInfo.pPoolSizes = &poolSize;

	descriptorPool = device.createDescriptorPool(createInfo);

	// Allocate the set in the pool

	vk::DescriptorSetAllocateInfo allocateInfo;
	allocateInfo.descriptorPool = descriptorPool;
	allocateInfo.descriptorSetCount = 1;
	allocateInfo.pSetLayouts = &descriptorSetLayout;

	descriptorSets = device.allocateDescriptorSets(allocateInfo);

	// Connect buffers with descriptor
	std::vector<vk::DescriptorBufferInfo> bufferInfos(buffers.size());
	std::vector<vk::WriteDescriptorSet> writeDescriptorSets(buffers.size());
	for (int i = 0; i < buffers.size(); i++) {
		bufferInfos[i].buffer = buffers[i];
		bufferInfos[i].offset = 0;
		bufferInfos[i].range = bufferSizes[i];

		writeDescriptorSets[i].dstSet = descriptorSets[0];
		writeDescriptorSets[i].dstBinding = i;
		writeDescriptorSets[i].descriptorCount = 1;
		writeDescriptorSets[i].descriptorType = vk::DescriptorType::eStorageBuffer;
		writeDescriptorSets[i].pBufferInfo = &bufferInfos[i];
	}

	device.updateDescriptorSets(writeDescriptorSets, nullptr);

}

void ComputeApp::createComputePipeline() {

	// Create Shader
	uint32_t filelength;
	std::vector<uint32_t> code = readShaderFile(filelength, "resources/comp.spv");
	vk::ShaderModuleCreateInfo shaderCreateInfo;
	shaderCreateInfo.codeSize = filelength;
	shaderCreateInfo.pCode = code.data();

	computeShaderModule = device.createShaderModule(shaderCreateInfo);

	// Specify compute shader stage
	vk::PipelineShaderStageCreateInfo shaderStageCreateInfo;
	shaderStageCreateInfo.stage = vk::ShaderStageFlagBits::eCompute;
	shaderStageCreateInfo.module = computeShaderModule;
	shaderStageCreateInfo.pName = "main";

	// Create PipelineLayout
	vk::PipelineLayoutCreateInfo pipelineLayoutCreateInfo;
	pipelineLayoutCreateInfo.setLayoutCount = 1;
	pipelineLayoutCreateInfo.pSetLayouts = &descriptorSetLayout;

	pipelineLayout = device.createPipelineLayout(pipelineLayoutCreateInfo);

	// Create Pipeline
	vk::ComputePipelineCreateInfo pipelineCreateInfo;
	pipelineCreateInfo.stage = shaderStageCreateInfo;
	pipelineCreateInfo.layout = pipelineLayout;

	pipeline = device.createComputePipeline(nullptr, pipelineCreateInfo);
}

void ComputeApp::createCommandeBuffer() {
	auto queueFamilyIndex = findQueueFamilies(physicalDevice);

	// Create Command Pool
	vk::CommandPoolCreateInfo poolCreateInfo;
	poolCreateInfo.queueFamilyIndex = static_cast<uint32_t>(queueFamilyIndex.computeFamily);
	commandPool = device.createCommandPool(poolCreateInfo);

	// Create Command buffer
	vk::CommandBufferAllocateInfo allocateInfo;
	allocateInfo.level = vk::CommandBufferLevel::ePrimary;
	allocateInfo.commandPool = commandPool;
	allocateInfo.commandBufferCount = 1;
	commandBuffer = device.allocateCommandBuffers(allocateInfo)[0];

	// Record some commands
	vk::CommandBufferBeginInfo beginInfo;
	beginInfo.flags = vk::CommandBufferUsageFlagBits::eOneTimeSubmit;
	commandBuffer.begin(beginInfo);

	commandBuffer.bindPipeline(vk::PipelineBindPoint::eCompute, pipeline);
	commandBuffer.bindDescriptorSets(vk::PipelineBindPoint::eCompute, pipelineLayout, /*first set*/0, descriptorSets, nullptr);

	// run the compute shader
	commandBuffer.dispatch(workgroupSize, 1, 1);

	commandBuffer.end();
}

bool ComputeApp::isDeviceSuitable(vk::PhysicalDevice device) {
	QueueFamilyIndices indices = findQueueFamilies(device);

	return indices.isComplete();
}

ComputeApp::QueueFamilyIndices ComputeApp::findQueueFamilies(vk::PhysicalDevice physDevice) {
	QueueFamilyIndices indices;

	auto queueFamilyProperties = physDevice.getQueueFamilyProperties();

	int i = 0;
	for (const auto& queueFamily : queueFamilyProperties) {
		if (queueFamily.queueCount > 0 && queueFamily.queueFlags & vk::QueueFlagBits::eCompute) {
			indices.computeFamily = i;
		}

		if (indices.isComplete())
			break;

		++i;
	}

	return indices;
}

uint32_t ComputeApp::findMemoryType(uint32_t memoryTypeBits, vk::MemoryPropertyFlags properties) {
	auto memoryProperties = physicalDevice.getMemoryProperties();

	for (uint32_t i = 0; i < memoryProperties.memoryTypeCount; i++) {
		if (
			(memoryTypeBits & (i << 1)) &&
			((memoryProperties.memoryTypes[i].propertyFlags & properties) == properties)
		)
			return i;
	}
	return -1;
}
