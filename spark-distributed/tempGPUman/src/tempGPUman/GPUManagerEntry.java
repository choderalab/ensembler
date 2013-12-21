package tempGPUman;

public class GPUManagerEntry {
	
	private gpuManager gpuman;
	
	public GPUManagerEntry(int numGPU, int workersPerGPU)
	{
		gpuman = new gpuManager(numGPU, workersPerGPU);
	}
	
	public gpuManager getgpuMan()
	{
		return gpuman;
	}

}
