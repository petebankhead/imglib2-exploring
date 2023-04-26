/**
 * QuPath script to explore getting an ImgLib2 image backed by
 * a QuPath ImageServer.
 *
 * (Current implementation assumes that all pixels can be cached)
 *
 * Also, dependencies need to be available (as specified in build.gradle).
 */



import groovy.transform.CompileStatic
import net.imglib2.RandomAccessibleInterval
import net.imglib2.img.Img
import net.imglib2.img.basictypeaccess.array.ByteArray
import net.imglib2.img.basictypeaccess.array.FloatArray
import net.imglib2.img.basictypeaccess.array.IntArray
import net.imglib2.img.basictypeaccess.array.ShortArray
import net.imglib2.img.cell.Cell
import net.imglib2.img.cell.CellGrid
import net.imglib2.img.cell.LazyCellImg
import net.imglib2.img.display.imagej.ImageJFunctions
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory
import net.imglib2.realtransform.RealViews
import net.imglib2.realtransform.Scale2D
import net.imglib2.type.NativeType
import net.imglib2.type.numeric.ARGBType
import net.imglib2.type.numeric.integer.*
import net.imglib2.type.numeric.real.DoubleType
import net.imglib2.type.numeric.real.FloatType
import net.imglib2.view.IntervalView
import net.imglib2.view.Views
import org.scijava.util.DoubleArray
import qupath.lib.images.servers.ImageServer
import qupath.lib.images.servers.ImageServers
import qupath.lib.images.servers.PixelType
import qupath.lib.images.servers.ServerTools
import qupath.lib.images.servers.TileRequest

import java.awt.image.BufferedImage
import java.awt.image.Raster
import java.util.concurrent.ConcurrentHashMap

def path = "/path/to/image.tif"

def server = ImageServers.buildServer(path)
def factory = wrapServer(server)
//def img = factory.createForLevel(server.nResolutions()-2)
def img = factory.createForMaxDimension(1034)
ImageJFunctions.show(img)


QuPathServerImgBuilder wrapServer(ImageServer<BufferedImage> server) {
    return new QuPathServerImgBuilder(
            server,
            getNativeType(server))
}

def <T> NativeType<T> getNativeType(ImageServer<?> server) {
    if (server.isRGB())
        return new ARGBType()
    switch (server.getPixelType()) {
        case PixelType.INT8:
            return new ByteType();
        case PixelType.UINT8:
            return new UnsignedByteType();
        case PixelType.INT16:
            return new ShortType();
        case PixelType.UINT16:
            return new UnsignedShortType();
        case PixelType.INT32:
            return new IntType();
        case PixelType.UINT32:
            return new UnsignedIntType();
        case PixelType.FLOAT32:
            return new FloatType();
        case PixelType.FLOAT64:
            return new DoubleType();
        throw new UnsupportedOperationException("Unsupported type " + server.getPixelType())
    }
}


class QuPathServerImgBuilder<T extends NativeType<T>> {

    private ImageServer<BufferedImage> server;
    private T nativeType;

    // Need to cache the cells or performance is bad
    private Map<Integer, Map<Long, Cell>> caches = new ConcurrentHashMap<>()

    QuPathServerImgBuilder(ImageServer<BufferedImage> server, T nativeType) {
        this.server = server;
        this.nativeType = nativeType;
    }

    private Map<Long, Cell> getCellCache(Integer level) {
        return caches.computeIfAbsent(level, ConcurrentHashMap::new)
    }

    Img<T> createForLevel(int level) {
        int nChannels = nativeType instanceof ARGBType ? 1 : server.nChannels()
        long[] dims = [server.getMetadata().getLevel(level).getWidth(),
                      server.getMetadata().getLevel(level).getHeight(),
                       nChannels] as long[]
        int[] cellDims = [server.getMetadata().getPreferredTileWidth(),
                          server.getMetadata().getPreferredTileHeight(),
                          nChannels] as int[]
        var tiles = server
                .getTileRequestManager()
                .getTileRequestsForLevel(level) as List
        var cache = getCellCache(level);
        return new LazyCellImg(
                new CellGrid(dims, cellDims),
                nativeType,
                (long v) -> cache.computeIfAbsent(v, i -> getCell(server, tiles[(int)i]))
        )
    }

    RandomAccessibleInterval<T> createForDownsample(double downsample) {
        int level = ServerTools.getPreferredResolutionLevel(server, downsample);
        def img = createForLevel(level);
        double downsampleForLevel = server.getDownsampleForResolution(level);
        double scaleFactor = downsampleForLevel / downsample;
        return scale(img, scaleFactor);
    }

    RandomAccessibleInterval<T> createForMaxDimension(int maxLength) {
        // TODO: Improve this - it can give off-by-one results
        double downsample = Math.max(server.getWidth(), server.getHeight()) / maxLength
        return createForDownsample(downsample);
    }

    private Cell getCell(ImageServer<BufferedImage> server, TileRequest tile) {
        def img = server.readRegion(tile.getRegionRequest())
        int nChannels = nativeType instanceof ARGBType ? 1 : server.nChannels()
        int[] dims = [img.getWidth(), img.getHeight(), nChannels] as int[]
        long[] min = [tile.getTileX(), tile.getTileY(), 0] as long[]
        def data
        if (nativeType instanceof FloatType)
            data = getFloats(img.getRaster())
        else if (nativeType instanceof DoubleType)
            data = getDoubles(img.getRaster())
        else if (nativeType instanceof ARGBType)
            data = getARGB(img)
        else if (nativeType instanceof UnsignedIntType || nativeType instanceof IntType)
            data = getIntegers(img.getRaster())
        else if (nativeType instanceof UnsignedByteType || nativeType instanceof ByteType)
            data = getBytes(img.getRaster())
        else if (nativeType instanceof UnsignedShortType || nativeType instanceof ShortType)
            data = getShorts(img.getRaster())
        else
            throw new IllegalArgumentException("Unsupported type " + nativeType)
        return new Cell<>(
                dims,
                min,
                data
        )
    }

    @CompileStatic
    private static FloatArray getFloats(Raster raster) {
        int planeSize = getPlaneSize(raster);
        int nBands = raster.getNumBands();
        float[] array = new float[planeSize * nBands];
        float[] source = null;
        for (int band = 0; band < nBands; band++) {
            source = raster.getSamples(0, 0, raster.getWidth(), raster.getHeight(), band, source);
            System.arraycopy(source, 0, array, band * planeSize, planeSize)
        }
        return new FloatArray(array);
    }

    @CompileStatic
    private static DoubleArray getDoubles(Raster raster) {
        int planeSize = getPlaneSize(raster);
        int nBands = raster.getNumBands();
        double[] array = new double[planeSize * nBands];
        double[] source = null;
        for (int band = 0; band < nBands; band++) {
            source = raster.getSamples(0, 0, raster.getWidth(), raster.getHeight(), band, source);
            System.arraycopy(source, 0, array, band * planeSize, planeSize)
        }
        return new DoubleArray(array);
    }

    @CompileStatic
    private static ByteArray getBytes(Raster raster) {
        int planeSize = getPlaneSize(raster);
        int nBands = raster.getNumBands();
        byte[] array = new byte[planeSize * nBands];
        int[] source = null;
        int ind = 0;
        for (int band = 0; band < nBands; band++) {
            source = raster.getSamples(0, 0, raster.getWidth(), raster.getHeight(), band, source);
            for (int val : source) {
                array[ind++] = (byte)val;
            }
        }
        return new ByteArray(array);
    }

    @CompileStatic
    private static ShortArray getShorts(Raster raster) {
        int planeSize = getPlaneSize(raster);
        int nBands = raster.getNumBands();
        short[] array = new short[planeSize * nBands];
        int[] source = null;
        int ind = 0;
        for (int band = 0; band < nBands; band++) {
            source = raster.getSamples(0, 0, raster.getWidth(), raster.getHeight(), band, source);
            for (int val : source) {
                array[ind++] = (short)val;
            }
        }
        return new ShortArray(array);
    }

    @CompileStatic
    private static IntArray getIntegers(Raster raster) {
        int planeSize = getPlaneSize(raster);
        int nBands = raster.getNumBands();
        int[] array = new byte[planeSize * nBands];
        int[] source = null;
        for (int band = 0; band < nBands; band++) {
            source = raster.getSamples(0, 0, raster.getWidth(), raster.getHeight(), band, source);
            System.arraycopy(source, 0, array, band * planeSize, planeSize);
        }
        return new IntArray(array);
    }

    @CompileStatic
    private static IntArray getARGB(BufferedImage img) {
        int[] argb = img.getRGB(0, 0, img.getWidth(), img.getHeight(), null, 0, img.getWidth());
        return new IntArray(argb);
    }

    private static int getPlaneSize(Raster raster) {
        return raster.getWidth() * raster.getHeight();
    }

    /**
     * Scale an image
     * @param img
     * @param scaleFactor
     * @return
     */
    static <T> RandomAccessibleInterval<T> scale(RandomAccessibleInterval<T> img, double scaleFactor) {
        if (scaleFactor == 1.0)
            return img;
        else if (scaleFactor <= 0 || !Double.isFinite(scaleFactor))
            throw new IllegalArgumentException("Invalid scale factor " + scaleFactor + " - must be a positive real number")
        long[] min = img.minAsLongArray()
        long[] max = img.maxAsLongArray()
        max[0] = Math.round(max[0] * scaleFactor)
        max[1] = Math.round(max[1] * scaleFactor)
        def scale = new Scale2D(scaleFactor, scaleFactor)
        return Views.interval(
                Views.raster(
                        RealViews.affine(
                                Views.interpolate(
                                        Views.extendMirrorDouble(img),
                                        new NLinearInterpolatorFactory<>()),
                                scale),
                ), min, max
        )
    }

}